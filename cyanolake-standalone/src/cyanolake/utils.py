"""Standalone utils implementation."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

from . import settings as _m_settings


# Reference cell 12, lines 1-2.
def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


# Reference cell 12, lines 5-10.
def human_size(size: int | float) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024


# Reference cell 12, lines 13-18.
def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


# Reference cell 12, lines 21-27.
def write_json(path: Path, payload: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    temporary.replace(path)
    return path


# Reference cell 12, lines 30-64.
def run_logged(
    command: list[str], log_path: Path, timeout: int | None = None, cwd: Path | None = None
) -> subprocess.CompletedProcess:
    """Run a command, write a complete log, and raise with a useful tail."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    result = subprocess.run(
        [str(item) for item in command],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=_m_settings.OCSSW_SUBPROCESS_ENV.copy(),
    )
    elapsed = time.time() - started
    header = f"started_utc={utc_now()}\nelapsed_seconds={elapsed:.3f}\nreturn_code={result.returncode}\ncommand={shlex.join([str(item) for item in command])}\n\n"
    log_path.write_text(
        header
        + "===== STDOUT =====\n"
        + (result.stdout or "")
        + "\n===== STDERR =====\n"
        + (result.stderr or ""),
        encoding="utf-8",
    )
    if result.returncode != 0:
        combined = ((result.stdout or "") + "\n" + (result.stderr or "")).strip()
        raise RuntimeError(
            f"Command failed with return code {result.returncode}. Log: {log_path}\n\n{combined[-4000:]}"
        )
    return result
