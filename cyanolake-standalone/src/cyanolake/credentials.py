"""Explicit credential prompts; secrets never enter configuration or provenance."""

from __future__ import annotations

import getpass
import netrc
import os
import tempfile
from pathlib import Path


def cdse_credentials() -> tuple[str, str]:
    username = os.environ.get("CDSE_USERNAME") or input("CDSE username: ").strip()
    password = os.environ.get("CDSE_PASSWORD") or getpass.getpass("CDSE password (hidden): ")
    if not username or not password:
        raise ValueError("CDSE username and password are required")
    return username, password


def configure_earthdata() -> Path:
    username = input("NASA Earthdata username: ").strip()
    password = getpass.getpass("NASA Earthdata password (hidden): ")
    if not username or not password or any(c.isspace() for c in username + password):
        raise ValueError("Use nonempty Earthdata credentials without whitespace")
    path = Path.home() / ".netrc"
    existing = path.read_text() if path.exists() else ""
    if path.is_symlink():
        raise ValueError("Refusing to replace a symbolic-link .netrc")
    if existing:
        # Preserve unrelated machines and macros exactly; append only new entries.
        parsed = netrc.netrc(str(path))
        conflicts = [
            h
            for h in ("urs.earthdata.nasa.gov", "oceandata.sci.gsfc.nasa.gov")
            if parsed.authenticators(h)
        ]
        if conflicts:
            raise ValueError(
                "Earthdata entries already exist; edit those entries manually to update them"
            )

    def quote(token):
        return '"' + token.replace("\\", "\\\\").replace('"', '\\"') + '"'

    content = existing.rstrip() + "\n"
    for host in ("urs.earthdata.nasa.gov", "oceandata.sci.gsfc.nasa.gov"):
        content += f"machine {host}\n  login {quote(username)}\n  password {quote(password)}\n"
    fd, tmp = tempfile.mkstemp(prefix=".netrc-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(content)
        netrc.netrc(tmp)
        os.replace(tmp, path)
    finally:
        Path(tmp).unlink(missing_ok=True)
    return path
