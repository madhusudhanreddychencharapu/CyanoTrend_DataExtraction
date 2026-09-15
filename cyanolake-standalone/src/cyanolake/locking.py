"""One process owns a workspace during writes; the OS releases locks on exit."""

import fcntl
from contextlib import contextmanager

from . import settings


@contextmanager
def workspace_lock():
    path = settings.PERSISTENT_ROOT / ".workspace.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "Another process is using this workspace; use a separate workspace or wait"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
