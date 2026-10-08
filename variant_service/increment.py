"""Incremental callset updates: sample replacement, deterministic re-merge, diff.

Pure logic on top of the existing parsing/normalization/merging pipeline:

* an increment only carries *new* or *corrected* samples; each carried sample
  wholesale replaces that sample's previous calls (a correction is never
  stacked on top of old conclusions as extra evidence);
* entries/conflicts at sites no increment touches are reused byte-for-byte, so
  untouched entries never disappear or reorder;
* the version id is a hash of the full materialized callset, so the same
  baseline plus the same increment always yields the same version and summary.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json

from .merging import _conflict_key, _sort_key, merge_callset
from .models import (
    CallsetResult,
    ChangeSummary,
    Conflict,
    ConflictChange,
    EntryChange,
    EntryChangeType,
    Genotype,
    Rejection,
    SampleChange,
    SampleChangeType,
    VariantRecord,
    VersionedCallset,
)


class BaselineNotFound(KeyError):
    """Raised when an increment references an unknown/evicted base version."""

    code = "BASELINE_NOT_FOUND"


def _loc_tuple(loc) -> tuple:
    return (loc.block, loc.line_no, loc.raw)


def _conflict_tuple(c: Conflict) -> tuple:
    return (
        c.conflict_type.value,
        c.variant.render(),
        c.sample,
        c.detail,
        [_loc_tuple(l) for l in c.locations],
    )


def canonical_payload(result: CallsetResult) -> bytes:
    """Deterministic JSON serialization of the materialized callset.

    Only callset *content* participates: ``content_hash`` (derived from the
    request text and shape) deliberately does not, so equivalent callsets
    reached via different submission paths share one version id.
    """
    doc = {
        "entries": [
            {
                "variant": e.variant.render(),
                "samples": list(e.samples),
                "genotype_by_sample": {s: g.value for s, g in e.genotype_by_sample.items()},
                "support": [_loc_tuple(l) for l in e.support],
                "conflicts": [_conflict_tuple(c) for c in e.conflicts],
            }
            for e in result.entries
        ],
        "rejected": [
            {
                "reason": r.reason.value,
                "detail": r.detail,
                "sample": r.sample,
                "location": _loc_tuple(r.location),
            }
            for r in result.rejected
        ],
        "conflicts": [_conflict_tuple(c) for c in result.conflicts],
        "accepted_count": result.accepted_count,
    }
    return json.dumps(doc, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def canonical_version_id(result: CallsetResult) -> str:
    """Hash the full materialized callset content into a stable version id."""
    return hashlib.sha256(canonical_payload(result)).hexdigest()


def _partition_rejections(
    rejected: list[Rejection],
) -> tuple[dict[str, list[Rejection]], list[Rejection]]:
    """Split rejects into sample-attributable and batch-level (sample=None)."""
    by_sample: dict[str, list[Rejection]] = {}
    batch: list[Rejection] = []
    for r in rejected:
        if r.sample is None:
            batch.append(r)
        else:
            by_sample.setdefault(r.sample, []).append(r)
    return by_sample, batch


def _ordered_samples(accepted, rejected_by_sample: dict, sample_order: tuple[str, ...]) -> list[str]:
    found = {r.sample for r in accepted} | set(rejected_by_sample)
    ordered = [s for s in sample_order if s in found]
    known = set(ordered)
    ordered.extend(sorted(found - known))
    return ordered


def _materialize(
    content_hash: str,
    accepted_by_sample: dict[str, tuple[VariantRecord, ...]],
    rejected_by_sample: dict[str, tuple[Rejection, ...]],
    batch_rejected: dict[frozenset[str], tuple[Rejection, ...]],
) -> CallsetResult:
    """Flatten per-sample state into one accepted/rejected stream and re-merge.

    Batch-level rejects are emitted in deterministic batch-key order, then
    sample-attributable rejects, preserving the existing merge output contract.
    """
    accepted: list[VariantRecord] = []
    for sample in sorted(accepted_by_sample):
        accepted.extend(accepted_by_sample[sample])
    rejected: list[Rejection] = []
    for key in sorted(batch_rejected, key=lambda k: tuple(sorted(k))):
        rejected.extend(batch_rejected[key])
    for sample in sorted(rejected_by_sample):
        rejected.extend(rejected_by_sample[sample])
    return merge_callset(content_hash, accepted, rejected)


def _entry_sample_changes(before, after) -> tuple[SampleChange, ...]:
    """Per-sample carrier/genotype differences for one retained variant."""
    changes: list[SampleChange] = []
    b_gts = before.genotype_by_sample
    a_gts = after.genotype_by_sample
    for sample in sorted(set(b_gts) | set(a_gts)):
        g_before = b_gts.get(sample)
        g_after = a_gts.get(sample)
        bv = g_before.value if g_before is not None else None
        av = g_after.value if g_after is not None else None
        was_carrier = g_before is not None and g_before in (Genotype.HET, Genotype.HOM)
        now_carrier = g_after is not None and g_after in (Genotype.HET, Genotype.HOM)
        if bv == av:
            continue
        if not was_carrier and now_carrier:
            changes.append(SampleChange(SampleChangeType.CARRIER_ADDED, sample, bv, av))
        elif was_carrier and not now_carrier:
            changes.append(SampleChange(SampleChangeType.CARRIER_REMOVED, sample, bv, av))
        else:
            changes.append(SampleChange(SampleChangeType.GENOTYPE_CHANGED, sample, bv, av))
    return tuple(changes)


def build_change_summary(
    base: VersionedCallset | None,
    current: VersionedCallset,
    delta_samples: tuple[str, ...],
) -> ChangeSummary:
    """Diff a baseline (or None for an initial build) against ``current``."""
    before_entries = {e.variant: e for e in base.result.entries} if base else {}
    after_entries = {e.variant: e for e in current.result.entries}

    added: list[str] = []
    removed: list[str] = []
    modified: list[EntryChange] = []
    for key in after_entries:
        rendered = key.render()
        if key not in before_entries:
            added.append(rendered)
        else:
            changes = _entry_sample_changes(before_entries[key], after_entries[key])
            before = before_entries[key]
            after = after_entries[key]
            if changes or before.samples != after.samples:
                modified.append(EntryChange(
                    change_type=EntryChangeType.MODIFIED,
                    variant=rendered,
                    samples_before=before.samples,
                    samples_after=after.samples,
                    sample_changes=changes,
                ))
    for key in before_entries:
        if key not in after_entries:
            removed.append(key.render())

    def conflict_ids(result):
        out = {}
        for c in result.conflicts:
            out[_conflict_key(c)] = c
        return out

    before_c = conflict_ids(base.result) if base else {}
    after_c = conflict_ids(current.result)
    appeared = [
        ConflictChange(True, c.conflict_type, c.variant.render(), c.sample, c.detail)
        for k, c in after_c.items() if k not in before_c
    ]
    resolved = [
        ConflictChange(False, c.conflict_type, c.variant.render(), c.sample, c.detail)
        for k, c in before_c.items() if k not in after_c
    ]
    appeared.sort(key=lambda x: (x.variant, x.conflict_type.value, x.sample or "", x.detail))
    resolved.sort(key=lambda x: (x.variant, x.conflict_type.value, x.sample or "", x.detail))
    added.sort()
    removed.sort()
    modified.sort(key=lambda x: x.variant)

    base_samples = set(base.sample_order) if base else set()
    delta_set = frozenset(delta_samples)
    new_samples = tuple(sorted(delta_set - base_samples))
    corrected = tuple(s for s in delta_samples if s in base_samples)

    return ChangeSummary(
        baseline_version=base.version_id if base else None,
        version=current.version_id,
        submitted_samples=tuple(dict.fromkeys(delta_samples)),
        corrected_samples=corrected,
        new_samples=new_samples,
        entries_added=tuple(added),
        entries_removed=tuple(removed),
        entry_changes=tuple(modified),
        conflicts_appeared=tuple(appeared),
        conflicts_resolved=tuple(resolved),
    )


def build_initial_version(
    content_hash: str,
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    sample_order: tuple[str, ...],
) -> VersionedCallset:
    """Build the first immutable version from a full submission."""
    rej_by_sample, batch = _partition_rejections(rejected)

    accepted_by_sample: dict[str, tuple[VariantRecord, ...]] = {}
    grouped: dict[str, list[VariantRecord]] = {}
    for rec in accepted:
        grouped.setdefault(rec.sample, []).append(rec)
    for sample, recs in grouped.items():
        accepted_by_sample[sample] = tuple(recs)

    rejected_by_sample = {s: tuple(rs) for s, rs in rej_by_sample.items()}
    ordered = _ordered_samples(accepted, rejected_by_sample, sample_order)
    batch_rejected = {frozenset(ordered): tuple(batch)} if batch else {}

    result = _materialize(content_hash, accepted_by_sample, rejected_by_sample, batch_rejected)
    version_id = canonical_version_id(result)
    version = VersionedCallset(
        version_id=version_id,
        parent_version_id=None,
        result=result,
        accepted_by_sample=accepted_by_sample,
        rejected_by_sample=rejected_by_sample,
        batch_rejected=batch_rejected,
        sample_order=tuple(ordered),
    )
    summary = build_change_summary(None, version, tuple(ordered))
    return dataclasses.replace(version, summary=summary)


def _site_of(rec: VariantRecord) -> tuple[str, int]:
    return (rec.key.chrom, rec.key.pos)


def _touched_scope(
    old_accepted: dict[str, tuple[VariantRecord, ...]],
    new_accepted: dict[str, tuple[VariantRecord, ...]],
    delta: frozenset[str],
) -> tuple[set[tuple[str, int]], set]:
    """Sites/variant keys potentially affected by replacing the delta samples."""
    sites: set[tuple[str, int]] = set()
    for sample in delta:
        for rec in old_accepted.get(sample, ()):
            sites.add(_site_of(rec))
        for rec in new_accepted.get(sample, ()):
            sites.add(_site_of(rec))
    # every variant anchored at a touched site participates (another sample's
    # call at that site may form SAMPLE_CALL_CONFLICT/REF_CONFLICT with the
    # corrected sample's new conclusion)
    keys: set = set()
    for pool in (old_accepted, new_accepted):
        for recs in pool.values():
            for rec in recs:
                if _site_of(rec) in sites:
                    keys.add(rec.key)
    return sites, keys


def _scoped_merge(
    content_hash: str,
    accepted_by_sample: dict[str, tuple[VariantRecord, ...]],
    touched_keys: set,
) -> tuple[list, list]:
    """Run the shared merger over records restricted to touched variant keys."""
    accepted: list[VariantRecord] = []
    for sample in sorted(accepted_by_sample):
        for rec in accepted_by_sample[sample]:
            if rec.key in touched_keys:
                accepted.append(rec)
    scoped = merge_callset(content_hash, accepted, [])
    return scoped.entries, scoped.conflicts


def apply_delta(
    base: VersionedCallset,
    content_hash: str,
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    delta_samples: tuple[str, ...],
) -> VersionedCallset:
    """Apply an increment (whole-sample replacement) onto ``base``."""
    delta = frozenset(delta_samples)

    # 1) Wholesale-replace per-sample accepted calls and sample-attributable
    #    rejects: old conclusions of a corrected sample are removed first, so
    #    they are never stacked with the new submission as extra evidence.
    accepted_by_sample: dict[str, tuple[VariantRecord, ...]] = {
        s: recs for s, recs in base.accepted_by_sample.items() if s not in delta
    }
    rejected_by_sample: dict[str, tuple[Rejection, ...]] = {
        s: rs for s, rs in base.rejected_by_sample.items() if s not in delta
    }
    new_grouped: dict[str, list[VariantRecord]] = {}
    for rec in accepted:
        new_grouped.setdefault(rec.sample, []).append(rec)
    for sample, recs in new_grouped.items():
        accepted_by_sample[sample] = tuple(recs)
    new_rej_by_sample: dict[str, list[Rejection]] = {}
    batch_rej: list[Rejection] = []
    for r in rejected:
        if r.sample is None:
            batch_rej.append(r)
        else:
            new_rej_by_sample.setdefault(r.sample, []).append(r)
    for sample in delta:
        new_rs = tuple(new_rej_by_sample.get(sample, ()))
        if new_rs:
            rejected_by_sample[sample] = new_rs

    # 2) Batch-level structural rejects: this increment batch replaces the
    #    rejects previously owned by exactly the same declared sample set;
    #    rejects of every other historical batch stay untouched.
    batch_rejected: dict[frozenset[str], tuple[Rejection, ...]] = dict(base.batch_rejected)
    if batch_rej:
        batch_rejected[delta] = tuple(batch_rej)
    elif delta in batch_rejected:
        del batch_rejected[delta]

    # 3) Stable sample order: historical first-seen order preserved, brand-new
    #    samples appended in deterministic header order.
    present = set(accepted_by_sample) | set(rejected_by_sample)
    known = set(base.sample_order)
    order = [s for s in base.sample_order if s in present]
    for sample in delta_samples:
        if sample in present and sample not in known and sample not in order:
            order.append(sample)
    for sample in sorted(present - set(order)):
        order.append(sample)

    # 4) Targeted re-merge of affected variant keys; untouched entries are
    #    reused byte-for-byte and never reordered or dropped.
    _, touched_keys = _touched_scope(base.accepted_by_sample, accepted_by_sample, delta)
    rebuilt_entries, rebuilt_conflicts = _scoped_merge(
        content_hash, accepted_by_sample, touched_keys
    )
    rebuilt_by_key = {e.variant: e for e in rebuilt_entries}
    final_entries: list = []
    for entry in base.result.entries:
        if entry.variant in touched_keys:
            if entry.variant in rebuilt_by_key:
                final_entries.append(rebuilt_by_key.pop(entry.variant))
        else:
            final_entries.append(entry)
    for entry in sorted(rebuilt_by_key.values(), key=lambda e: _sort_key(e.variant)):
        final_entries.append(entry)
    # Re-establish the merger's global deterministic order. Untouched entries
    # keep object identity; new/rebuilt entries land at their sorted position.
    final_entries.sort(key=lambda e: _sort_key(e.variant))

    # 5) Conflicts: old ones attached to touched variants are replaced by the
    #    rebuilt set; conflicts at untouched sites remain exactly as before.
    rebuilt_ids = {_conflict_key(c) for c in rebuilt_conflicts}
    final_conflicts = [
        c for c in base.result.conflicts
        if c.variant not in touched_keys and _conflict_key(c) not in rebuilt_ids
    ]
    final_conflicts.extend(rebuilt_conflicts)
    final_conflicts.sort(key=_conflict_key)

    # 6) Assemble the materialized rejected surface and the new result.
    flat_rejected: list[Rejection] = []
    for key in sorted(batch_rejected, key=lambda k: tuple(sorted(k))):
        flat_rejected.extend(batch_rejected[key])
    for sample in sorted(rejected_by_sample):
        flat_rejected.extend(rejected_by_sample[sample])
    accepted_count = sum(len(recs) for recs in accepted_by_sample.values())
    result = CallsetResult(
        content_hash=content_hash,
        entries=final_entries,
        rejected=flat_rejected,
        conflicts=final_conflicts,
        accepted_count=accepted_count,
    )

    version = VersionedCallset(
        version_id=canonical_version_id(result),
        parent_version_id=base.version_id,
        result=result,
        accepted_by_sample=accepted_by_sample,
        rejected_by_sample=rejected_by_sample,
        batch_rejected=batch_rejected,
        sample_order=tuple(order),
        delta_samples=delta,
    )
    summary = build_change_summary(base, version, delta_samples)
    return dataclasses.replace(version, summary=summary)
