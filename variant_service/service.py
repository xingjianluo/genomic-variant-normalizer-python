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

import dataclasses
import hashlib
import logging
import threading
from collections import OrderedDict

from .increment import build_increment, build_root, restate_against
from .merging import merge_callset, validate_and_normalize
from .models import (
    CallsetResult,
    CallsetSnapshot,
    IncrementErrorType,
    RawRecord,
    Rejection,
)
from .parsing import parse_multi_sample, parse_single_sample

logger = logging.getLogger(__name__)

_CACHE_MAX = 128
_VERSION_MAX = 128


class BaselineError(Exception):
    """An incremental commit cannot be applied to the requested baseline."""

    def __init__(
        self,
        error_type: IncrementErrorType,
        detail: str,
        *,
        current_version: str | None = None,
        conflicting_samples: tuple[str, ...] = (),
    ) -> None:
        super().__init__(detail)
        self.error_type = error_type
        self.detail = detail
        self.current_version = current_version
        self.conflicting_samples = conflicting_samples


class CallsetService:
    """Stateless-per-request processing with a shared content-addressed cache."""

    def __init__(self, cache_max: int = _CACHE_MAX, version_max: int = _VERSION_MAX) -> None:
        self._cache: OrderedDict[str, CallsetResult] = OrderedDict()
        self._cache_max = cache_max
        # version_id -> immutable snapshot; plus per-lineage bookkeeping.
        self._versions: OrderedDict[str, CallsetSnapshot] = OrderedDict()
        self._version_max = version_max
        # root_version -> {"head": str, "touched": {sample: version_id}}
        self._lineage: dict[str, dict] = {}
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

    @staticmethod
    def _owners_for_result(raw: list[RawRecord], rejected: list[Rejection]):
        """Default rejection ownership for a stateless legacy submission."""
        line_samples: dict[int, tuple[str, ...]] = {}
        for rr in raw:
            line_samples[rr.location.line_no] = tuple(
                sorted(set(line_samples.get(rr.location.line_no, ())) | {rr.sample})
            )
        return [line_samples.get(r.location.line_no, ()) for r in rejected]

    def _register_result_as_root(self, result: CallsetResult,
                                 raw: list[RawRecord]) -> CallsetSnapshot:
        """Register a legacy submission's result as a retrievable root version.

        Legacy endpoints remain unchanged, but their result gains a stable
        ``version_id`` clients can later cite as an increment baseline. Rebuild
        is avoided by reconstructing the snapshot's record state directly from
        the merged result is lossy, so the snapshot is derived from raw inputs
        via the normal funnel and content-addressed (idempotent -> same object).
        """
        accepted, rejected = validate_and_normalize(raw)
        owners = self._owners_for_result(raw, rejected)
        snap = build_root(accepted, rejected, owners)
        self._store_snapshot(snap)
        return snap

    def _run_cached(self, content_hash: str, raw: list[RawRecord],
                    rejected: list[Rejection],
                    register_version: bool = False) -> CallsetResult:
        # Fast path under lock for cache hits; compute without holding the lock
        # so a heavy request never blocks cache access for other requests.
        with self._lock:
            cached = self._cache.get(content_hash)
            if cached is not None:
                self._cache.move_to_end(content_hash)
                logger.info("提交命中缓存 hash=%s 前 8 位，复用既有结论 entries=%d",
                            content_hash[:8], len(cached.entries))
                if register_version and cached.version_id is None:
                    snap = self._register_result_as_root(cached, raw)
                    cached = dataclasses.replace(cached, version_id=snap.version_id)
                    self._cache[content_hash] = cached
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
        if register_version:
            # Make the submission citable as an increment baseline. The root is
            # content-addressed independently of the submission content hash.
            snap = self._register_result_as_root(result, raw)
            result = dataclasses.replace(result, version_id=snap.version_id)
        return result

    def submit_single(self, text: str, sample: str) -> CallsetResult:
        sample = sample.strip()
        content_hash = self.hash_content(text, sample)
        cached = self.get_cached(content_hash)
        if cached is not None and cached.version_id is not None:
            return cached
        raw, rejected = parse_single_sample(text, sample)
        return self._run_cached(content_hash, raw, rejected, register_version=True)

    def submit_multi(self, text: str) -> CallsetResult:
        content_hash = self.hash_content(text)
        cached = self.get_cached(content_hash)
        if cached is not None and cached.version_id is not None:
            return cached
        raw, rejected, _names = parse_multi_sample(text)
        return self._run_cached(content_hash, raw, rejected, register_version=True)

    # ------------------------------------------------------------------
    # Versioned incremental updates
    # ------------------------------------------------------------------

    def get_version(self, version_id: str) -> CallsetSnapshot | None:
        """Return an immutable snapshot by version id, or None (read under lock).

        The lookup is O(1) and only takes the lock long enough to fetch from an
        OrderedDict, so a large in-flight rebuild never blocks version reads.
        """
        with self._lock:
            snap = self._versions.get(version_id)
            if snap is not None:
                self._versions.move_to_end(version_id)
            return snap

    def _normalize_payload(
        self, raw: list[RawRecord], rejected: list[Rejection]
    ) -> tuple[list, list[Rejection]]:
        """Run the shared semantic validation funnel for an increment payload."""
        accepted, norm_rejected = validate_and_normalize(raw)
        return accepted, rejected + norm_rejected

    def _store_snapshot_locked(self, snap: CallsetSnapshot,
                               commit_seq: int | None = None) -> int:
        """Register snapshot and advance lineage head. Caller holds the lock.

        Returns the monotonic lineage commit sequence assigned to this store.
        """
        if commit_seq is None:
            lineage = self._lineage.get(snap.root_version)
            commit_seq = (lineage["commit_seq"] + 1) if lineage else 0
        if snap.commit_seq != commit_seq:
            snap = dataclasses.replace(snap, commit_seq=commit_seq)
        if snap.version_id not in self._versions:
            self._versions[snap.version_id] = snap
            self._versions.move_to_end(snap.version_id)
            while len(self._versions) > self._version_max:
                victim_id, _victim = next(iter(self._versions.items()))
                # Never evict a lineage head: it is still a valid baseline.
                if (
                    snap.root_version in self._lineage
                    and self._lineage[snap.root_version]["head"] == victim_id
                ):
                    self._versions.move_to_end(victim_id)
                    break
                self._versions.popitem(last=False)
        lineage = self._lineage.setdefault(
            snap.root_version,
            {"head": snap.version_id, "touched": {}, "commit_seq": commit_seq - 1},
        )
        lineage["commit_seq"] = commit_seq
        lineage["head"] = snap.version_id
        for sample in snap.touched_samples:
            if commit_seq > lineage["touched"].get(sample, -1):
                lineage["touched"][sample] = commit_seq
        return commit_seq

    def _store_snapshot(self, snap: CallsetSnapshot) -> None:
        with self._lock:
            self._store_snapshot_locked(snap)

    def _commit_root(
        self, accepted: list, rejected: list[Rejection], owners: list[tuple[str, ...]]
    ) -> CallsetSnapshot:
        """Build, register and return a root snapshot (sequence 0).

        Content-addressed, so an identical root payload reuses the same snapshot
        object and is inherently idempotent.
        """
        snap = build_root(accepted, rejected, owners)
        self._store_snapshot(snap)
        stored = self.get_version(snap.version_id)
        logger.info(
            "建立根版本 version=%s 前 8 位: 样本 %d 个, 条目 %d, 拒绝 %d, 冲突 %d",
            snap.version_id[:8], len(snap.samples), len(snap.result.entries),
            len(snap.result.rejected), len(snap.result.conflicts),
        )
        return stored if stored is not None else snap

    def _commit_increment(
        self,
        base_version: str,
        accepted: list,
        rejected: list[Rejection],
        touched_samples: list[str],
        owners: list[tuple[str, ...]],
    ) -> CallsetSnapshot:
        """Optimistically commit an increment on top of a lineage head.

        Concurrency protocol (deterministic, no lost samples):

        * the baseline and current lineage head are read under the lock, but the
          rebuild runs *without* holding the lock so big rebuilds never block
          version reads or other requests;
        * at store time the head is re-checked under the lock (CAS): if another
          commit landed in between and it changed a sample this commit also
          touches, the two updates cannot be merged deterministically and the
          request is rejected as ``BASELINE_SUPERSEDED`` with the current head
          and conflicting samples — the client retries against the head;
        * if the intervening commit only touched *disjoint* samples, the
          increment is rebuilt once on top of the new head (rebase) and the CAS
          is retried; content addressing guarantees both commit orders converge
          to the exact same final version.
        """
        touched = set(touched_samples)
        if not touched:
            # A payload carrying no named sample at all (e.g. only malformed
            # columns) cannot replace anything; the API guards empty text.
            raise BaselineError(
                IncrementErrorType.EMPTY_INCREMENT,
                "增量没有携带任何样本，无法应用",
            )

        with self._lock:
            baseline = self._versions.get(base_version)
            if baseline is None:
                raise BaselineError(
                    IncrementErrorType.BASELINE_NOT_FOUND,
                    f"未知或已淘汰的基线版本标识: {base_version}",
                )
            lineage = self._lineage.get(baseline.root_version)
            if lineage is None:
                raise BaselineError(
                    IncrementErrorType.BASELINE_NOT_FOUND,
                    f"版本 {base_version} 的血缘已不存在",
                )

        # Idempotency probe computed WITHOUT the lock so it never blocks reads
        # of existing versions: what would this increment yield if applied
        # straight onto the cited baseline?
        wanted = build_increment(baseline, accepted, rejected, sorted(touched), owners)
        if wanted.version_id == baseline.version_id:
            # No change relative to the cited baseline: genuine no-op.
            return wanted

        with self._lock:
            head_id = lineage["head"]
            head_snap = self._versions.get(head_id)
            known = self._versions.get(wanted.version_id)

        def _is_ancestor_or_head(candidate_id: str | None) -> bool:
            """True if the candidate state is already part of the head lineage.

            Head itself, or an earlier committed state (smaller commit_seq): a
            replay of either produces content the lineage already contains, so
            no newer correction is discarded.
            """
            if candidate_id is None or head_snap is None:
                return False
            if candidate_id == head_snap.version_id:
                return True
            candidate = self._versions.get(candidate_id)
            return candidate is not None and candidate.commit_seq <= lineage["commit_seq"]

        if known is not None and _is_ancestor_or_head(known.version_id):
            # Replaying an ancestor/head state: the resulting content is already
            # part of the lineage, so no newer correction is discarded.
            logger.info(
                "增量命中既有版本 version=%s 前 8 位 (幂等/祖先重放)，按基线 %s 前 8 位重述摘要",
                known.version_id[:8], base_version[:8],
            )
            return restate_against(known, baseline)

        with self._lock:
            overlapping_now = sorted(touched & {
                s for s, seq in lineage["touched"].items()
                if seq > baseline.commit_seq
            })

        # Overlapping samples changed by a newer commit, and the wanted state is
        # not already contained in the head: accepting would discard the newer
        # correction or rewind the lineage. Deterministic, distinguishable reject.
        if overlapping_now:
            raise BaselineError(
                IncrementErrorType.BASELINE_SUPERSEDED,
                "基线已被含重叠样本的更新取代，无法确定性合并；请基于当前版本重试",
                current_version=head_id,
                conflicting_samples=tuple(overlapping_now),
            )

        max_attempts = 32
        current_head_id = head_id
        snap: CallsetSnapshot | None = None
        for attempt in range(max_attempts):
            with self._lock:
                head_at_read = lineage["head"]
                # Re-check overlap against any commit that landed while the
                # idempotency probe was running lock-free.
                changed_now = {
                    s for s, seq in lineage["touched"].items()
                    if seq > baseline.commit_seq
                }
            overlap_now = sorted(touched & changed_now)
            if overlap_now:
                raise BaselineError(
                    IncrementErrorType.BASELINE_SUPERSEDED,
                    "提交期间基线被含重叠样本的更新抢占；请基于当前版本重试",
                    current_version=lineage["head"],
                    conflicting_samples=tuple(overlap_now),
                )

            head = self.get_version(head_at_read)
            if head is None:  # evicted while rebuilding -> retry against live head
                current_head_id = lineage["head"]
                continue
            snap = build_increment(head, accepted, rejected, sorted(touched), owners)

            with self._lock:
                if lineage["head"] == head_at_read:
                    # CAS won: no intervening commit. Assign the next monotonic
                    # commit sequence and persist the rebased snapshot.
                    next_seq = lineage["commit_seq"] + 1
                    snap = dataclasses.replace(snap, commit_seq=next_seq)
                    self._store_snapshot_locked(snap, commit_seq=next_seq)
                    stored = self._versions[snap.version_id]
                    break
                # Someone committed first; loop and rebase onto the new head.
                current_head_id = lineage["head"]
                logger.info(
                    "增量 CAS 竞争 base=%s: 头版本已前进到 %s 前 8 位，第 %d 次重放",
                    base_version[:8], current_head_id[:8], attempt + 1,
                )
        else:  # pragma: no cover - bounded retries exhausted under extreme load
            raise BaselineError(
                IncrementErrorType.BASELINE_SUPERSEDED,
                "头版本竞争过于激烈，多次重放仍未追上；请稍后基于当前版本重试",
                current_version=current_head_id,
            )

        logger.info(
            "增量提交 base=%s -> version=%s 前 8 位 (seq=%d): 新增条目 %d, 变化 %d, "
            "移除 %d, 新冲突 %d, 消除冲突 %d, identical=%s",
            base_version[:8], stored.version_id[:8], stored.sequence,
            len(stored.changes.added_entries), len(stored.changes.changed_entries),
            len(stored.changes.removed_entries), len(stored.changes.new_conflicts),
            len(stored.changes.resolved_conflicts), stored.changes.identical,
        )
        stored2 = self.get_version(stored.version_id)
        return stored2 if stored2 is not None else stored

    @staticmethod
    def _owners_for(
        names: list[str], raw: list[RawRecord], rejected: list[Rejection]
    ) -> list[tuple[str, ...]]:
        """Attribute structural rejections of a multi-sample payload.

        A header-level problem (duplicate names, missing sample columns) applies
        to every named sample. A data-row problem is attributed to the sample
        columns that row supplied when those raw records exist; otherwise the
        whole payload owns it.
        """
        all_names = tuple(sorted(names))
        samples_by_line: dict[int, tuple[str, ...]] = {}
        for rr in raw:
            samples_by_line.setdefault(rr.location.line_no, ())
            samples_by_line[rr.location.line_no] = tuple(
                sorted(set(samples_by_line[rr.location.line_no]) | {rr.sample})
            )
        owners: list[tuple[str, ...]] = []
        for r in rejected:
            owners.append(samples_by_line.get(r.location.line_no) or all_names)
        return owners

    def submit_root_single(self, text: str, sample: str) -> CallsetSnapshot:
        """Create a root lineage from a single-sample payload."""
        sample = sample.strip()
        raw, rejected = parse_single_sample(text, sample)
        accepted, rejected = self._normalize_payload(raw, rejected)
        return self._commit_root(accepted, rejected, [(sample,)] * len(rejected))

    def submit_root_multi(self, text: str) -> CallsetSnapshot:
        """Create a root lineage from a multi-sample payload."""
        raw, rejected, names = parse_multi_sample(text)
        accepted, rejected = self._normalize_payload(raw, rejected)
        return self._commit_root(accepted, rejected, self._owners_for(names, raw, rejected))

    def submit_increment_single(
        self, text: str, sample: str, base_version: str
    ) -> CallsetSnapshot:
        """Apply one single-sample payload as an increment to base_version.

        Empty (but non-whitespace-protected at the API) text means the sample
        retracts all its calls; the sample name itself is always the touched
        key, so an empty payload still replaces (clears) that sample.
        """
        sample = sample.strip()
        raw, rejected = parse_single_sample(text, sample)
        accepted, rejected = self._normalize_payload(raw, rejected)
        return self._commit_increment(
            base_version, accepted, rejected, [sample], [(sample,)] * len(rejected)
        )

    def submit_increment_multi(self, text: str, base_version: str) -> CallsetSnapshot:
        """Apply a multi-sample payload as an increment to base_version."""
        raw, rejected, names = parse_multi_sample(text)
        if not names:
            # No parseable sample columns: nothing to add or replace.
            raise BaselineError(
                IncrementErrorType.EMPTY_INCREMENT,
                "增量文本未包含任何样本列，无法应用",
            )
        accepted, rejected = self._normalize_payload(raw, rejected)
        return self._commit_increment(
            base_version, accepted, rejected, names,
            self._owners_for(names, raw, rejected),
        )
