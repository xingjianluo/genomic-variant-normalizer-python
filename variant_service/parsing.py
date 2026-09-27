"""Parse single-sample and multi-sample variant text into the shared model."""
from __future__ import annotations

from .models import RawRecord, Rejection


def parse_single_sample(text: str, sample: str) -> tuple[list[RawRecord], list[Rejection]]:
    raise NotImplementedError


def parse_multi_sample(text: str) -> tuple[list[RawRecord], list[Rejection]]:
    raise NotImplementedError
