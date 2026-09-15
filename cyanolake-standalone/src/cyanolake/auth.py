"""Standalone auth implementation."""

from __future__ import annotations

import time

import requests

from . import settings as _m_settings


# Reference cell 14, lines 9-36.
class CDSEAuth:
    def __init__(self) -> None:
        self.access_token: str | None = None
        self.refresh_token: str | None = None
        self.expires_at = 0.0

    def login(self, username: str, password: str) -> str:
        response = requests.post(
            _m_settings.CDSE_TOKEN_URL,
            data={
                "client_id": "cdse-public",
                "grant_type": "password",
                "username": username,
                "password": password,
            },
            timeout=60,
        )
        response.raise_for_status()
        payload = response.json()
        self.access_token = payload["access_token"]
        self.refresh_token = payload.get("refresh_token")
        self.expires_at = time.time() + int(payload.get("expires_in", 600)) - 60
        return self.access_token

    def token(self, username: str, password: str) -> str:
        if self.access_token and time.time() < self.expires_at:
            return self.access_token
        return self.login(username, password)


CDSE_AUTH = CDSEAuth()
