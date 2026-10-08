"""Incremental callset updates: per-sample replace, re-merge, versioning, diff.

Pure functions over immutable :class:`CallsetSnapshot` state; no I/O, no locks.
Persistence and concurrency control live in :mod:`variant_service.service`.

State model
-----------
A snapshot stores the *full* set of normalized records attributed to samples.
An increment never adds evidence to an existing sample: every touched sample's
records are first removed and then replaced wholesale by the new payload, so a
correction cannot stack two conclusions. Untouched samples are carried over
unchanged (same tuples, same ordering), and the callset is then rebuilt with
the exact same merge/conflict funnel used for a first submission.

Version identity
----------------
``version_id`` is ``sha256`` over a canonical serialization of the complete
record state. It is a *content* address (the parent pointer is lineage
metadata, not part of the identity): the same state always has the same
version, even when two disjoint increments land in different orders. Same
parent + same replacement payload therefore yields the same version and the
same :class:`ChangeSummary`, no matter how many times it is submitted.
"""
from __future__ import annotations

import hashlib

from .merging import _conflict_key, _loc_key, merge_callset
from .models import (
    CallsetEntry,
    CallsetResult,
    CallsetSnapshot,
    ChangeSummary,
    Conflict,
    EntryChange,
    Rejection,
    VariantKey,
    VariantRecord,
)

_STATE_VERSION_TAG = "callset-state-v1"


def _rejection_identity(r: Rejection) -> tuple:
    return (r.reason.value, *_loc_key(r.location), r.detail)


def _carriers(entry: CallsetEntry) -> set[str]:
    """Samples whose presence is visible on the entry (sorted carriers tuple)."""
    return set(entry.samples)


def _entry_map(result: CallsetResult) -> dict[VariantKey, CallsetEntry]:
    return {e.variant: e for e in result.entries}


def _diff_entries(old: CallsetResult, new: CallsetResult) -> tuple:
    """Return (added, removed, changed) entry key lists, deterministically."""
    old_map, new_map = _entry_map(old), _entry_map(new)
    added = sorted(new_map.keys() - old_map.keys(), key=lambda k: k.render())
    removed = sorted(old_map.keys() - new_map.keys(), key=lambda k: k.render())
    changed: list[EntryChange] = []
    for key in sorted(old_map.keys() & new_map.keys(), key=lambda k: k.render()):
        oe, ne = old_map[key], new_map[key]
        old_carriers, new_carriers = _carriers(oe), _carriers(ne)
        samples_added = tuple(sorted(new_carriers - old_carriers))
        samples_removed = tuple(sorted(old_carriers - new_carriers))
        genotype_changes: list[tuple[str, str, str]] = []
        for sample in sorted(set(oe.genotype_by_sample) | set(ne.genotype_by_sample)):
            old_g = oe.genotype_by_sample.get(sample)
            new_g = ne.genotype_by_sample.get(sample)
            if old_g is not new_g and (old_g is None or new_g is None or old_g != new_g):
                genotype_changes.append(
                    (sample, old_g.value if old_g else "./.", new_g.value if new_g else "./.")
                )
        if samples_added or samples_removed or genotype_changes:
            changed.append(
                EntryChange(
                    variant=key,
                    kind="CHANGED",
                    samples_added=samples_added,
                    samples_removed=samples_removed,
                    genotype_changes=tuple(genotype_changes),
                )
            )
    return added, removed, changed


def _all_carrier_samples(result: CallsetResult) -> set[str]:
    """Samples with at least one ALT carrier presence across entries."""
    out: set[str] = set()
    for e in result.entries:
        out.update(e.samples)
    return out


def _diff_conflicts(old: CallsetResult, new: CallsetResult) -> tuple:
    old_by = {_conflict_key(c): c for c in old.conflicts}
    new_by = {_conflict_key(c): c for c in new.conflicts}
    new_keys = sorted(new_by.keys() - old_by.keys())
    resolved_keys = sorted(old_by.keys() - new_by.keys())
    return (
        tuple(new_by[k] for k in new_keys),
        tuple(old_by[k] for k in resolved_keys),
    )


def _diff_rejections(
    old: CallsetResult,
    old_owners: tuple[tuple[str, ...], ...],
    new: CallsetResult,
    new_owners: tuple[tuple[str, ...], ...],
) -> tuple:
    old_by = {_rejection_identity(r): r for r in old.rejected}
    new_by = {_rejection_identity(r): r for r in new.rejected}
    new_keys = {_rejection_identity(r) for r in new.rejected}
    old_keys = {_rejection_identity(r) for r in old.rejected}
    return (
        tuple(sorted((new_by[k] for k in new_keys - old_keys), key=_rejection_identity)),
        tuple(sorted((old_by[k] for k in old_keys - new_keys), key=_rejection_identity)),
    )


def _diff(
    old: CallsetResult,
    old_owners: tuple[tuple[str, ...], ...],
    new: CallsetResult,
    new_owners: tuple[tuple[str, ...], ...],
) -> ChangeSummary:
    added, removed, changed = _diff_entries(old, new)
    new_conflicts, resolved_conflicts = _diff_conflicts(old, new)
    new_rejections, resolved_rejections = _diff_rejections(
        old, old_owners, new, new_owners
    )
    old_carriers, new_carriers = _all_carrier_samples(old), _all_carrier_samples(new)
    samples_added = tuple(sorted(new_carriers - old_carriers))
    samples_removed = tuple(sorted(old_carriers - new_carriers))
    identical = not (
        added or removed or changed
        or samples_added or samples_removed
        or new_conflicts or resolved_conflicts
        or new_rejections or resolved_rejections
    )
    return ChangeSummary(
        added_entries=tuple(added),
        removed_entries=tuple(removed),
        changed_entries=tuple(changed),
        samples_added=samples_added,
        samples_removed=samples_removed,
        new_conflicts=new_conflicts,
        resolved_conflicts=resolved_conflicts,
        new_rejections=new_rejections,
        resolved_rejections=resolved_rejections,
        identical=identical,
    )


def _canonical_accepted(rec: VariantRecord) -> str:
    return "\t".join(
        (
            rec.sample,
            rec.key.chrom,
            str(rec.key.pos),
            rec.key.ref,
            rec.key.alt,
            rec.genotype.value,
            rec.location.block,
            str(rec.location.line_no),
            rec.location.raw,
        )
    )


def _canonical_rejected(r: Rejection, owners: tuple[str, ...]) -> str:
    return "\t".join(
        (
            ",".join(owners),
            r.reason.value,
            r.location.block,
            str(r.location.line_no),
            r.location.raw,
            r.detail,
        )
    )


def _state_hash(
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    owners: list[tuple[str, ...]],
) -> str:
    lines = [_STATE_VERSION_TAG]
    lines.extend(_canonical_accepted(r) for r in sorted(accepted, key=_canonical_accepted))
    lines.extend(
        _canonical_rejected(r, o)
        for r, o in sorted(zip(rejected, owners), key=lambda ro: _canonical_rejected(*ro))
    )
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def _rebuild(
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    content_hash: str,
) -> CallsetResult:
    return merge_callset(content_hash, list(accepted), list(rejected))


def _sorted_parallel(rejected, owners):
    pairs = sorted(zip(rejected, owners), key=lambda ro: _loc_key(ro[0].location))
    return (
        tuple(r for r, _ in pairs),
        tuple(o for _, o in pairs),
    )


def build_root(
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    owners: list[tuple[str, ...]] | None = None,
) -> CallsetSnapshot:
    """Build the root snapshot (sequence 0) for a first full submission."""
    if owners is None:
        owners = [tuple(sorted({r.location.block.removeprefix("sample:")}))
                  if r.location.block.startswith("sample:") else ()
                  for r in rejected]
    rejected_s, owners_s = _sorted_parallel(rejected, owners)
    accepted_s = tuple(sorted(accepted, key=_canonical_accepted))
    version_id = _state_hash(list(accepted_s), list(rejected_s), list(owners_s))
    result = _rebuild(list(accepted_s), list(rejected_s), version_id)
    samples = tuple(sorted({r.sample for r in accepted_s} | {s for o in owners_s for s in o}))
    return CallsetSnapshot(
        version_id=version_id,
        parent_version=None,
        root_version=version_id,
        sequence=0,
        commit_seq=0,
        samples=samples,
        touched_samples=samples,
        accepted=accepted_s,
        rejected=rejected_s,
        rejection_owners=owners_s,
        result=result,
        changes=ChangeSummary(identical=False),
    )


def _owner_index(rejected, owners):
    return {_rejection_identity(r): tuple(sorted(o)) for r, o in zip(rejected, owners)}


def restate_against(
    snap: CallsetSnapshot, baseline: CallsetSnapshot
) -> CallsetSnapshot:
    """Return the same state as ``snap`` with a change summary vs ``baseline``.

    The version id and callset content are untouched (content addressing), but
    the reported :class:`ChangeSummary` is computed against the baseline the
    client actually referenced. Used when a state already exists from another
    lineage path but is re-requested against a different parent.
    """
    changes = _diff(
        baseline.result, baseline.rejection_owners,
        snap.result, snap.rejection_owners,
    )
    return CallsetSnapshot(
        version_id=snap.version_id,
        parent_version=baseline.version_id,
        root_version=snap.root_version,
        sequence=snap.sequence,
        commit_seq=snap.commit_seq,
        samples=snap.samples,
        touched_samples=snap.touched_samples,
        accepted=snap.accepted,
        rejected=snap.rejected,
        rejection_owners=snap.rejection_owners,
        result=snap.result,
        changes=changes,
    )


def build_increment(
    baseline: CallsetSnapshot,
    accepted: list[VariantRecord],
    rejected: list[Rejection],
    touched_samples: list[str],
    owners: list[tuple[str, ...]] | None = None,
) -> CallsetSnapshot:
    """Replace touched samples wholesale on top of baseline and re-merge.

    Records (accepted and rejected) previously attributed to any touched sample
    are removed first; the payload then supplies that sample's entire new
    conclusion. Evidence is never stacked across versions.
    """
    touched = set(touched_samples)
    if owners is None:
        owners = [tuple(sorted(touched))] * len(rejected)

    kept_accepted = [r for r in baseline.accepted if r.sample not in touched]
    old_owner_by_id = _owner_index(baseline.rejected, baseline.rejection_owners)
    kept_rej: list[Rejection] = []
    kept_owners: list[tuple[str, ...]] = []
    for r in baseline.rejected:
        o = old_owner_by_id.get(_rejection_identity(r), ())
        if not (set(o) & touched):
            kept_rej.append(r)
            kept_owners.append(o)

    all_accepted = kept_accepted + list(accepted)
    all_rejected = kept_rej + list(rejected)
    all_owners = kept_owners + [tuple(sorted(o)) for o in owners]

    rejected_s, owners_s = _sorted_parallel(all_rejected, all_owners)
    accepted_s = tuple(sorted(all_accepted, key=_canonical_accepted))
    version_id = _state_hash(list(accepted_s), list(rejected_s), list(owners_s))
    result = _rebuild(list(accepted_s), list(rejected_s), version_id)
    samples = tuple(sorted({r.sample for r in accepted_s} | {s for o in owners_s for s in o}))

    # An increment whose full state equals the parent's state keeps the parent
    # version: resubmitting the same payload on the same baseline is idempotent.
    if version_id == baseline.version_id:
        return CallsetSnapshot(
            version_id=baseline.version_id,
            parent_version=baseline.version_id,
            root_version=baseline.root_version,
            sequence=baseline.sequence,
            commit_seq=baseline.commit_seq,
            samples=baseline.samples,
            touched_samples=tuple(sorted(touched)),
            accepted=baseline.accepted,
            rejected=baseline.rejected,
            rejection_owners=baseline.rejection_owners,
            result=baseline.result,
            changes=ChangeSummary(identical=True),
        )

    changes = _diff(baseline.result, baseline.rejection_owners, result, owners_s)
    return CallsetSnapshot(
        version_id=version_id,
        parent_version=baseline.version_id,
        root_version=baseline.root_version,
        sequence=baseline.sequence + 1,
        commit_seq=baseline.commit_seq + 1,
        samples=samples,
        touched_samples=tuple(sorted(touched)),
        accepted=accepted_s,
        rejected=rejected_s,
        rejection_owners=owners_s,
        result=result,
        changes=changes,
    )
