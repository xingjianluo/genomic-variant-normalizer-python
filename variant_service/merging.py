"""Validate/normalize raw records and merge samples into one callset."""
from __future__ import annotations

from .models import CallsetResult, RawRecord, Rejection, VariantRecord


def validate_and_normalize(
    raw: list[RawRecord],
) -> tuple[list[VariantRecord], list[Rejection]]:
    raise NotImplementedError


def merge_callset(
    content_hash: str,
    accepted: list[VariantRecord],
    rejected: list[Rejection],
) -> CallsetResult:
    raise NotImplementedError
