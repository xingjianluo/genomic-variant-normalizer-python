"""Application service: submission processing, caching, bounded concurrency.

Concurrency model
-----------------
* every submission is processed from local data only; results live inside the
  returned :class:`CallsetResult` and no per-sample global state is kept, so two
  concurrent requests over overlapping sample sets cannot contaminate each
  other;
* identical submissions (same text + sample identity) hash to the same key and
  reuse the previously computed result. The cache is a bounded, thread-safe
  LRU guarded by an :class:`~threading.RLock`;
* the parsing/normalization work itself is CPU-bound and O(n log n) in the
  number of records, so a large request cannot wedge the process state held by
  other requests — it only spends ordinary interpreter time.
"""
from __future__ import annotations

import hashlib
import logging
import threading
from collections import OrderedDict

from .merging import merge_callset, validate_and_normalize
from .models import CallsetResult, RawRecord, Rejection
from .parsing import parse_multi_sample, parse_single_sample

logger = logging.getLogger(__name__)

_CACHE_MAX = 128


class CallsetService:
    """Stateless-per-request processing with a shared content-addressed cache."""

    def __init__(self, cache_max: int = _CACHE_MAX) -> None:
        self._cache: OrderedDict[str, CallsetResult] = OrderedDict()
        self._cache_max = cache_max
        self._lock = threading.RLock()

    @staticmethod
    def hash_content(text: str, sample: str | None = None) -> str:
        payload = (sample or "") + "\x00" + text
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def get_cached(self, content_hash: str) -> CallsetResult | None:
        with self._lock:
            hit = self._cache.get(content_hash)
            if hit is not None:
                self._cache.move_to_end(content_hash)
            return hit

    def _store(self, content_hash: str, result: CallsetResult) -> None:
        with self._lock:
            self._cache[content_hash] = result
            self._cache.move_to_end(content_hash)
            while len(self._cache) > self._cache_max:
                self._cache.popitem(last=False)

    def cache_info(self) -> dict[str, int]:
        with self._lock:
            return {"entries": len(self._cache), "max": self._cache_max}

    def clear_cache(self) -> None:
        with self._lock:
            self._cache.clear()

    @staticmethod
    def _process(raw: list[RawRecord], rejected: list[Rejection], content_hash: str) -> CallsetResult:
        accepted, norm_rejected = validate_and_normalize(raw)
        all_rejected = rejected + norm_rejected
        return merge_callset(content_hash, accepted, all_rejected)

    def _run_cached(self, content_hash: str, raw: list[RawRecord],
                    rejected: list[Rejection]) -> CallsetResult:
        # Fast path under lock for cache hits; compute without holding the lock
        # so a heavy request never blocks cache access for other requests.
        with self._lock:
            cached = self._cache.get(content_hash)
            if cached is not None:
                self._cache.move_to_end(content_hash)
                logger.info("提交命中缓存 hash=%s 前 8 位，复用既有结论 entries=%d",
                            content_hash[:8], len(cached.entries))
                return cached

        logger.info("新提交 hash=%s 前 8 位: 原始记录 %d 条, 结构性拒绝 %d 条",
                    content_hash[:8], len(raw), len(rejected))
        result = self._process(raw, rejected, content_hash)
        logger.info(
            "判定依据 hash=%s: 接受 %d, 聚合条目 %d, 拒绝 %d, 冲突 %d",
            content_hash[:8], result.accepted_count, len(result.entries),
            len(result.rejected), len(result.conflicts),
        )
        self._store(content_hash, result)
        return result

    def submit_single(self, text: str, sample: str) -> CallsetResult:
        sample = sample.strip()
        content_hash = self.hash_content(text, sample)
        cached = self.get_cached(content_hash)
        if cached is not None:
            return cached
        raw, rejected = parse_single_sample(text, sample)
        return self._run_cached(content_hash, raw, rejected)

    def submit_multi(self, text: str) -> CallsetResult:
        content_hash = self.hash_content(text)
        cached = self.get_cached(content_hash)
        if cached is not None:
            return cached
        raw, rejected = parse_multi_sample(text)
        return self._run_cached(content_hash, raw, rejected)
