"""Standalone staging implementation."""

from __future__ import annotations

import concurrent.futures as cf
import shutil
import time
import zipfile
from pathlib import Path
from typing import Any

import boto3
import requests
from tqdm.auto import tqdm

from . import auth as _m_auth
from . import settings as _m_settings
from . import state as _m_state
from . import utils as _m_utils


# Reference cell 15, lines 1-22.
def safe_extract_zip(zip_path: Path, destination: Path) -> Path:
    """Extract a ZIP after rejecting path traversal and validate one .SEN3 product."""
    zip_path = Path(zip_path)
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        bad_member = archive.testzip()
        if bad_member:
            raise RuntimeError(f"ZIP CRC validation failed at {bad_member}")
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if target != destination and destination not in target.parents:
                raise RuntimeError(f"Unsafe ZIP member: {member.filename}")
        archive.extractall(destination)
    products = sorted(
        {
            path.parent
            for path in destination.rglob("xfdumanifest.xml")
            if path.parent.name.endswith(".SEN3")
        }
    )
    if len(products) != 1:
        raise RuntimeError(f"Expected one extracted .SEN3 product, found {len(products)}")
    return products[0]


# Reference cell 15, lines 25-33.
def validate_sen3(path: Path) -> Path:
    path = Path(path).expanduser().resolve()
    manifest = path / "xfdumanifest.xml"
    if not path.is_dir() or not manifest.is_file():
        raise FileNotFoundError(f"Not a valid OLCI .SEN3 directory: {path}")
    radiance_files = list(path.glob("Oa*_radiance.nc"))
    if len(radiance_files) < 20:
        raise RuntimeError(
            f"Incomplete OLCI product: found only {len(radiance_files)} radiance files"
        )
    return path


# Reference cell 15, lines 36-101.
def download_cdse_scene(scene: dict[str, Any], username: str, password: str) -> tuple[Path, Path]:
    """Resume a CDSE ZIP download, preserve it, and stage its .SEN3 product."""
    product_name = str(scene["name"]).removesuffix(".SEN3")
    zip_path = _m_settings.RAW_L1_DIR / f"{product_name}.zip"
    partial = zip_path.with_suffix(".zip.part")
    stage_root = _m_settings.STAGING_DIR / product_name
    existing_products = list(stage_root.rglob("xfdumanifest.xml")) if stage_root.exists() else []
    if zip_path.is_file():
        with zipfile.ZipFile(zip_path) as archive:
            bad_member = archive.testzip()
            if bad_member:
                raise RuntimeError(f"Existing ZIP failed CRC at {bad_member}: {zip_path}")
        if existing_products:
            return (zip_path, validate_sen3(existing_products[0].parent))
        if stage_root.exists():
            shutil.rmtree(stage_root)
        return (zip_path, validate_sen3(safe_extract_zip(zip_path, stage_root)))
    token = _m_auth.CDSE_AUTH.token(username, password)
    headers = {"Authorization": f"Bearer {token}"}
    mode = "wb"
    if partial.exists() and partial.stat().st_size:
        headers["Range"] = f"bytes={partial.stat().st_size}-"
        mode = "ab"
    url = f"{_m_settings.CDSE_DOWNLOAD_ROOT}/Products({scene['id']})/$value"
    with requests.get(url, headers=headers, stream=True, timeout=(60, 900)) as response:
        if response.status_code == 200 and mode == "ab":
            mode = "wb"
        response.raise_for_status()
        initial = partial.stat().st_size if mode == "ab" and partial.exists() else 0
        remaining = int(response.headers.get("content-length") or 0)
        with (
            partial.open(mode) as stream,
            tqdm(
                total=initial + remaining if remaining else None,
                initial=initial,
                unit="B",
                unit_scale=True,
                desc=product_name[:32],
            ) as progress,
        ):
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    stream.write(chunk)
                    progress.update(len(chunk))
    partial.replace(zip_path)
    with zipfile.ZipFile(zip_path) as archive:
        bad_member = archive.testzip()
        if bad_member:
            raise RuntimeError(f"Downloaded ZIP failed CRC at {bad_member}")
    if stage_root.exists():
        shutil.rmtree(stage_root)
    sen3_path = safe_extract_zip(zip_path, stage_root)
    validate_sen3(sen3_path)
    archive_manifest = {
        "created_utc": _m_utils.utc_now(),
        "scene": scene,
        "raw_zip": str(zip_path),
        "raw_zip_size_bytes": zip_path.stat().st_size,
        "raw_zip_sha256": _m_utils.sha256_file(zip_path) if _m_settings.COMPUTE_SHA256 else None,
        "staged_sen3": str(sen3_path),
    }
    _m_utils.write_json(
        _m_settings.EXPORT_DIR / product_name / "l1_archive_manifest.json", archive_manifest
    )
    return (zip_path, sen3_path)


# Reference cell 15, lines 104-120.
def stage_existing_input(path: str | Path) -> tuple[Path | None, Path]:
    """Use an existing .SEN3 directory or stage an existing CDSE ZIP."""
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        return (None, validate_sen3(path))
    if path.suffix.lower() == ".zip" and path.is_file():
        product_name = path.stem
        stage_root = _m_settings.STAGING_DIR / product_name
        manifests = list(stage_root.rglob("xfdumanifest.xml")) if stage_root.exists() else []
        if manifests:
            sen3 = validate_sen3(manifests[0].parent)
        else:
            if stage_root.exists():
                shutil.rmtree(stage_root)
            sen3 = safe_extract_zip(path, stage_root)
        return (path, sen3)
    raise FileNotFoundError("Existing input must be a CDSE ZIP or an extracted .SEN3 directory")


# Reference cell 40, lines 9-36.
def temporary_s3_client(username: str, password: str):
    cached = _m_state._S3_CLIENT_CACHE.get(str(username))
    if cached and time.time() - cached[0] < 45 * 60:
        return cached[1]
    token = _m_auth.CDSE_AUTH.token(username, password)
    response = requests.post(
        _m_settings.CDSE_S3_KEYS_URL, headers={"Authorization": f"Bearer {token}"}, timeout=90
    )
    response.raise_for_status()
    data = response.json()
    access = data.get("access_id") or data.get("access_key") or data.get("access")
    secret = data.get("secret") or data.get("secret_key")
    if not access or not secret:
        raise RuntimeError("CDSE S3 credential service did not return access_id/secret")
    time.sleep(5)
    client = boto3.client(
        "s3",
        endpoint_url=_m_settings.CDSE_S3_ENDPOINT,
        aws_access_key_id=access,
        aws_secret_access_key=secret,
        region_name="default",
    )
    _m_state._S3_CLIENT_CACHE[str(username)] = (time.time(), client, access)
    return client


# Reference cell 40, lines 39-52.
def release_temporary_s3_credentials(username: str, password: str) -> None:
    """Best-effort deletion of a temporary S3 key created by this runtime."""
    cached = _m_state._S3_CLIENT_CACHE.pop(str(username), None)
    if not cached or len(cached) < 3:
        return
    access_id = cached[2]
    try:
        token = _m_auth.CDSE_AUTH.token(username, password)
        requests.delete(
            f"{_m_settings.CDSE_S3_KEYS_URL}/access_id/{access_id}",
            headers={"Authorization": f"Bearer {token}"},
            timeout=60,
        )
    except Exception as exc:
        print(f"Warning: temporary S3 key cleanup failed: {exc}")


# Reference cell 40, lines 55-62.
def _split_s3_path(s3_path: str) -> tuple[str, str]:
    text = str(s3_path or "").strip()
    if not text:
        raise ValueError("Catalogue record has no S3Path")
    parts = text.lstrip("/").split("/", 1)
    if len(parts) != 2:
        raise ValueError(f"Unexpected S3Path: {s3_path}")
    return (parts[0], parts[1].rstrip("/") + "/")


# Reference cell 40, lines 65-102.
def download_scene_s3_direct(
    scene: dict[str, Any], username: str, password: str, workers: int = 6
) -> Path:
    """Download the SEN3 object tree directly. Existing same-size files are reused."""
    client = temporary_s3_client(username, password)
    bucket, prefix = _split_s3_path(scene.get("s3_path"))
    target = _m_settings.STAGING_DIR / f"{scene['name']}.SEN3"
    target.mkdir(parents=True, exist_ok=True)
    paginator = client.get_paginator("list_objects_v2")
    objects = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        objects.extend(page.get("Contents", []))
    objects = [obj for obj in objects if not str(obj["Key"]).endswith("/")]
    if not objects:
        raise FileNotFoundError(f"No S3 objects found under {scene.get('s3_path')}")

    def fetch(obj):
        key = str(obj["Key"])
        relative = key[len(prefix) :]
        destination = (target / relative).resolve()
        if not destination.is_relative_to(target.resolve()):
            raise ValueError("S3 object path escapes its staging directory")
        destination.parent.mkdir(parents=True, exist_ok=True)
        expected = int(obj.get("Size") or 0)
        if destination.exists() and destination.stat().st_size == expected:
            return (expected, False)
        client.download_file(bucket, key, str(destination))
        return (expected, True)

    downloaded = 0
    with cf.ThreadPoolExecutor(max_workers=max(1, int(workers))) as executor:
        futures = [executor.submit(fetch, obj) for obj in objects]
        for future in tqdm(
            cf.as_completed(futures), total=len(futures), desc=f"S3 {scene['name'][:24]}"
        ):
            size, changed = future.result()
            if changed:
                downloaded += size
    validate_sen3(target)
    print(f"S3 direct stage ready: {target} (downloaded {_m_utils.human_size(downloaded)})")
    return target


# Reference cell 40, lines 105-120.
def stage_scene_efficient(
    scene: dict[str, Any], username: str, password: str, method: str = None
) -> tuple[Path | None, Path, str]:
    if method is None:
        method = _m_settings.DOWNLOAD_METHOD
    if scene.get("local_input"):
        raw, sen3 = stage_existing_input(scene["local_input"])
        return raw, sen3, "local_input"
    method = str(method).lower()
    if method not in {"auto", "s3", "zip"}:
        raise ValueError("download method must be auto, s3, or zip")
    if method in {"auto", "s3"}:
        try:
            return (None, download_scene_s3_direct(scene, username, password), "s3_direct")
        except Exception as exc:
            if method == "s3":
                raise
            print(f"S3 direct staging failed; using ZIP fallback: {exc}")
    raw_zip, sen3 = download_cdse_scene(scene, username, password)
    return (raw_zip, sen3, "odata_zip_fallback")
