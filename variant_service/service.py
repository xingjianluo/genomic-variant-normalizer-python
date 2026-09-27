"""Application service: submission processing, caching, bounded concurrency."""
from __future__ import annotations

import hashlib

from .models import CallsetResult
from .parsing import parse_multi_sample, parse_single_sample
from .merging import merge_callset, validate_and_normalize

_CACHE_MAX = 128


class CallsetService:
    def __init__(self, cache_max: int = _CACHE_MAX) -> None:
        self._cache: dict[str, CallsetResult] = {}
        self._cache_max = cache_max

    def submit_single(self, text: str, sample: str) -> CallsetResult:
        raise NotImplementedError

    def submit_multi(self, text: str) -> CallsetResult:
        raise NotImplementedError

    def _process(self, content_hash: str, raw, rejected) -> CallsetResult:
        accepted, norm_rejected = validate_and_normalize(raw)
        rejected = rejected + norm_rejected
        return merge_callset(content_hash, accepted, rejected)

    @staticmethod
    def hash_content(text: str, sample: str | None = None) -> str:
        payload = (sample or "") + "\x00" + text
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()
