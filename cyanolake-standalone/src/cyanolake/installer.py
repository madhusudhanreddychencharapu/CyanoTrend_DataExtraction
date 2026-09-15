"""Standalone installer implementation."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import urllib.request

from . import settings as _m_settings


# Reference cell 5, lines 24-28.
def operational_tag_key(tag: str) -> tuple[int, int]:
    m = re.fullmatch("V(\\d{4})\\.(\\d+)", str(tag).strip())
    if not m:
        raise ValueError(tag)
    return (int(m.group(1)), int(m.group(2)))


# Reference cell 5, lines 31-33.
def _tail(text: str, n: int = 80) -> str:
    lines = (text or "").splitlines()
    return "\n".join(lines[-n:])


# Reference cell 5, lines 36-43.
def _valid_ocssw_install() -> bool:
    required = (
        _m_settings.OCSSWROOT / "bin" / "l2gen",
        _m_settings.OCSSWROOT / "share" / "olci",
        _m_settings.OCSSWROOT / "share" / "olci" / "s3a",
        _m_settings.OCSSWROOT / "share" / "olci" / "s3b",
    )
    return all((p.exists() for p in required))


# Reference cell 5, lines 46-187.
def install_ocssw_for_olci() -> str:
    if _valid_ocssw_install():
        print(f"Reusing complete OCSSW installation at {_m_settings.OCSSWROOT}")
        return "existing installation"
    if not _m_settings.RUN_OCSSW_INSTALL:
        raise RuntimeError("RUN_OCSSW_INSTALL=False but a complete OCSSW installation is absent.")
    if _m_settings.OCSSWROOT.exists():
        print("An incomplete/partial OCSSW installation already exists.")
        raise RuntimeError(
            "Configured OCSSW root is incomplete; choose a new empty root or repair it manually"
        )
    _m_settings.OCSSW_INSTALLER_DIR.mkdir(parents=True, exist_ok=True)
    downloads = {
        "install_ocssw": "https://oceandata.sci.gsfc.nasa.gov/manifest/install_ocssw",
        "manifest.py": "https://oceandata.sci.gsfc.nasa.gov/manifest/manifest.py",
    }
    for filename, url in downloads.items():
        destination = _m_settings.OCSSW_INSTALLER_DIR / filename
        print(f"Downloading {filename} ...")
        urllib.request.urlretrieve(url, destination)
    installer = _m_settings.OCSSW_INSTALLER_DIR / "install_ocssw"
    installer.chmod(493)
    tags_proc = subprocess.run(
        [__import__("sys").executable, str(installer), "--list_tags"],
        cwd=_m_settings.OCSSW_INSTALLER_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    if tags_proc.returncode != 0:
        raise RuntimeError(
            "Could not query NASA OCSSW tags.\n\nSTDOUT:\n"
            + _tail(tags_proc.stdout)
            + "\n\nSTDERR:\n"
            + _tail(tags_proc.stderr)
        )
    operational_tags = sorted(
        {
            line.strip()
            for line in tags_proc.stdout.splitlines()
            if re.fullmatch("V\\d{4}\\.\\d+", line.strip())
        },
        key=operational_tag_key,
        reverse=True,
    )
    if not operational_tags:
        raise RuntimeError("NASA installer returned no operational OCSSW V tags.")
    candidates = []
    if _m_settings.PREFERRED_OCSSW_TAG:
        if _m_settings.PREFERRED_OCSSW_TAG in operational_tags:
            candidates.append(_m_settings.PREFERRED_OCSSW_TAG)
        else:
            print(
                f"Preferred tag {_m_settings.PREFERRED_OCSSW_TAG} is not listed; falling back to current operational tags."
            )
    for tag in operational_tags[: _m_settings.FALLBACK_TAG_COUNT]:
        if tag not in candidates:
            candidates.append(tag)
    print("OCSSW install candidates:", ", ".join(candidates))
    failures = []
    env = os.environ.copy()
    env["OCSSWROOT"] = str(_m_settings.OCSSWROOT)
    for tag in candidates:
        if _m_settings.OCSSWROOT.exists() and any(_m_settings.OCSSWROOT.iterdir()):
            import datetime

            stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
            _m_settings.OCSSWROOT.rename(
                _m_settings.OCSSWROOT.with_name(_m_settings.OCSSWROOT.name + ".failed-" + stamp)
            )
        _m_settings.OCSSWROOT.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            __import__("sys").executable,
            str(installer),
            "--install_dir",
            str(_m_settings.OCSSWROOT),
            "--tag",
            tag,
            "--seadas",
            "--olcis3a",
            "--olcis3b",
            "--verbose",
        ]
        print(f"\nInstalling OCSSW {tag} for Sentinel-3A/B OLCI ...")
        print("Command:", " ".join(cmd))
        proc = subprocess.run(
            cmd,
            cwd=_m_settings.OCSSW_INSTALLER_DIR,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        log_path = _m_settings.OCSSW_INSTALLER_DIR / f"install_{tag}.log"
        log_path.write_text(
            "COMMAND:\n"
            + " ".join(cmd)
            + "\n\nSTDOUT:\n"
            + (proc.stdout or "")
            + "\n\nSTDERR:\n"
            + (proc.stderr or ""),
            encoding="utf-8",
        )
        if proc.returncode == 0 and _valid_ocssw_install():
            os.environ["OCSSWROOT"] = str(_m_settings.OCSSWROOT)
            version = subprocess.run(
                [str(_m_settings.OCSSWROOT / "bin" / "l2gen"), "--version"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )
            print(f"OCSSW installation successful: {tag}")
            if version.stdout.strip() or version.stderr.strip():
                print("l2gen version:", (version.stdout or version.stderr).strip())
            print("Installer log:", log_path)
            return tag
        failure_summary = f"{tag}: exit={proc.returncode}\nSTDOUT tail:\n{_tail(proc.stdout)}\nSTDERR tail:\n{_tail(proc.stderr)}\nFull log: {log_path}"
        failures.append(failure_summary)
        print("\nInstallation attempt failed:")
        print(failure_summary)
    raise RuntimeError(
        "All OCSSW installation attempts failed. The detailed NASA installer output is shown above and saved under the configured OCSSW installer directory.\n\n"
        + "\n\n".join(failures)
    )
