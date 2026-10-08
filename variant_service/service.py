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

from .increment import (
    BaselineNotFound,
    apply_delta,
    build_initial_version,
)
from .merging import validate_and_normalize
from .models import CallsetResult, RawRecord, Rejection, VersionedCallset
from .parsing import declared_samples, parse_multi_sample, parse_single_sample

logger = logging.getLogger(__name__)

_CACHE_MAX = 128
_VERSION_MAX = 256


class EmptyIncrement(ValueError):
    """Raised when an increment carries no samples and no structural rejects."""

    code = "EMPTY_INCREMENT"


class CallsetService:
    """Stateless-per-request processing with a shared content-addressed cache."""

    def __init__(self, cache_max: int = _CACHE_MAX, version_max: int = _VERSION_MAX) -> None:
        self._cache: OrderedDict[str, CallsetResult] = OrderedDict()
        self._cache_max = cache_max
        self._lock = threading.RLock()
        # Immutable, content-addressed callset versions referenced by callers
        # via the returned ``version_id``. Versions never mutate, so reads of a
        # baseline never block on another request's heavy recomputation.
        self._versions: OrderedDict[str, VersionedCallset] = OrderedDict()
        self._version_max = version_max
        # request hash (baseline + payload) -> derived version id, for idempotency
        self._delta_index: dict[str, str] = {}
        # initial-submission content hash -> version id (result object identity)
        self._result_version: dict[str, str] = {}

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
                old_hash, _ = self._cache.popitem(last=False)
                self._result_version.pop(old_hash, None)

    def cache_info(self) -> dict[str, int]:
        with self._lock:
            return {"entries": len(self._cache), "max": self._cache_max}

    def clear_cache(self) -> None:
        with self._lock:
            self._cache.clear()
            self._versions.clear()
            self._delta_index.clear()
            self._result_version.clear()

    # ------------------------------------------------------------------
    # Version store
    # ------------------------------------------------------------------
    def get_version(self, version_id: str) -> VersionedCallset | None:
        """Fetch an immutable baseline version (short lock, never runs work)."""
        with self._lock:
            version = self._versions.get(version_id)
            if version is not None:
                self._versions.move_to_end(version_id)
            return version

    def _store_version(self, version: VersionedCallset) -> VersionedCallset:
        with self._lock:
            existing = self._versions.get(version.version_id)
            if existing is not None:
                # same content => same version; reuse the stored object so
                # repeated identical submissions converge to identity
                self._versions.move_to_end(version.version_id)
                return existing
            self._versions[version.version_id] = version
            self._versions.move_to_end(version.version_id)
            while len(self._versions) > self._version_max:
                old_id, _ = self._versions.popitem(last=False)
                self._delta_index = {
                    h: vid for h, vid in self._delta_index.items() if vid != old_id
                }
            return version

    def version_info(self) -> dict[str, int]:
        with self._lock:
            return {"versions": len(self._versions), "max": self._version_max}

    def _build_version(
        self,
        content_hash: str,
        raw: list[RawRecord],
        rejected: list[Rejection],
        sample_order: tuple[str, ...],
        base: VersionedCallset | None,
        delta_samples: tuple[str, ...],
    ) -> VersionedCallset:
        # Semantic validation/normalization and (re)merging run *outside* every
        # lock: CPU-bound O(n log n) work, local data only, so a large update
        # never wedges cache/version reads for other concurrent requests.
        accepted, norm_rejected = validate_and_normalize(raw)
        all_rejected = rejected + norm_rejected
        if base is None:
            version = build_initial_version(content_hash, accepted, all_rejected, sample_order)
        else:
            version = apply_delta(base, content_hash, accepted, all_rejected, delta_samples)
        logger.info(
            "判定依据 hash=%s: 接受 %d, 聚合条目 %d, 拒绝 %d, 冲突 %d, version=%s 前 10 位",
            content_hash[:8], version.result.accepted_count, len(version.result.entries),
            len(version.result.rejected), len(version.result.conflicts),
            version.version_id[:10],
        )
        stored = self._store_version(version)
        self._store(content_hash, stored.result)
        return stored

    def _run_version_cached(self, content_hash: str, raw: list[RawRecord],
                            rejected: list[Rejection],
                            sample_order: tuple[str, ...] = ()) -> VersionedCallset:
        # Fast path under lock: identical submissions reuse the same objects,
        # so repeated calls converge to identity (result is version is result).
        with self._lock:
            vid = self._result_version.get(content_hash)
        if vid is not None:
            hit = self.get_version(vid)
            if hit is not None:
                self.get_cached(content_hash)  # touch LRU order
                logger.info("提交命中版本缓存 hash=%s 前 8 位 version=%s，复用既有对象 entries=%d",
                            content_hash[:8], hit.version_id[:10], len(hit.result.entries))
                return hit

        logger.info("新提交 hash=%s 前 8 位: 原始记录 %d 条, 结构性拒绝 %d 条",
                    content_hash[:8], len(raw), len(rejected))
        version = self._build_version(content_hash, raw, rejected, sample_order, None, ())
        with self._lock:
            self._result_version[content_hash] = version.version_id
        return version

    def submit_single_version(self, text: str, sample: str) -> VersionedCallset:
        sample = sample.strip()
        content_hash = self.hash_content(text, sample)
        raw, rejected = parse_single_sample(text, sample)
        return self._run_version_cached(
            content_hash, raw, rejected, (sample,) if sample else ())

    def submit_multi_version(self, text: str) -> VersionedCallset:
        content_hash = self.hash_content(text)
        raw, rejected = parse_multi_sample(text)
        return self._run_version_cached(content_hash, raw, rejected, declared_samples(text))

    # Backward-compatible result-only entry points.
    def submit_single(self, text: str, sample: str) -> CallsetResult:
        return self.submit_single_version(text, sample).result

    def submit_multi(self, text: str) -> CallsetResult:
        return self.submit_multi_version(text).result

    # ------------------------------------------------------------------
    # Incremental updates
    # ------------------------------------------------------------------
    def _require_base(self, baseline_version: str) -> VersionedCallset:
        base = self.get_version(baseline_version)
        if base is None:
            raise BaselineNotFound(baseline_version)
        return base

    def _empty_delta(self, text: str, raw: list[RawRecord],
                     rejected: list[Rejection], samples: tuple[str, ...]) -> bool:
        if not text.strip():
            return True
        # header-only / rows carrying no sample column: nothing to replace
        return not samples and not raw and bool(rejected or text.strip())

    def increment_multi(self, text: str, baseline_version: str) -> VersionedCallset:
        """Apply a multi-sample (or single-sample shaped) increment."""
        base = self._require_base(baseline_version)
        raw, rejected = parse_multi_sample(text)
        samples = declared_samples(text)
        if self._empty_delta(text, raw, rejected, samples):
            raise EmptyIncrement("增量未携带任何样本（空文本或仅含无样本列的表头/数据）")
        content_hash = self.hash_content(text, "DELTA\x00" + baseline_version)
        return self._increment(base, content_hash, raw, rejected, samples)

    def increment_single(self, text: str, sample: str, baseline_version: str) -> VersionedCallset:
        """Apply a single-sample increment (one corrected/new sample)."""
        base = self._require_base(baseline_version)
        sample = sample.strip()
        raw, rejected = parse_single_sample(text, sample)
        if not text.strip() or not sample:
            raise EmptyIncrement("增量文本与样本名均不能为空")
        if not raw and not rejected:
            raise EmptyIncrement("增量未包含任何有效行或可报告的结构问题")
        content_hash = self.hash_content(text, "DELTA\x00" + baseline_version + "\x00" + sample)
        return self._increment(base, content_hash, raw, rejected, (sample,))

    def _increment(self, base: VersionedCallset, content_hash: str,
                   raw: list[RawRecord], rejected: list[Rejection],
                   samples: tuple[str, ...]) -> VersionedCallset:
        # Idempotent fast path: the same baseline + same payload maps to the
        # same request hash; reuse the previously derived version object.
        with self._lock:
            cached_vid = self._delta_index.get(content_hash)
        if cached_vid is not None:
            cached = self.get_version(cached_vid)
            if cached is not None:
                logger.info("增量命中幂等缓存 base=%s -> version=%s（复用既有版本对象）",
                            base.version_id[:8], cached.version_id[:10])
                return cached
        logger.info(
            "增量提交 base=%s 前 8 位 hash=%s 前 8 位: 样本 %s, 原始记录 %d, 结构拒绝 %d",
            base.version_id[:8], content_hash[:8], list(samples), len(raw), len(rejected),
        )
        version = self._build_version(content_hash, raw, rejected, samples, base, samples)
        with self._lock:
            self._delta_index[content_hash] = version.version_id
        summary = version.summary
        logger.info(
            "增量判定依据 version=%s: 新增条目 %d, 删除 %d, 修改 %d, 新冲突 %d, 消除冲突 %d",
            version.version_id[:10], len(summary.entries_added), len(summary.entries_removed),
            len(summary.entry_changes), len(summary.conflicts_appeared),
            len(summary.conflicts_resolved),
        )
        return version
