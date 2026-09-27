"""HTTP API exposing the offline service."""
from __future__ import annotations

from fastapi import APIRouter

from .service import CallsetService

router = APIRouter()
service = CallsetService()
