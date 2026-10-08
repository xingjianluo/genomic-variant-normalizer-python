"""Validate/normalize raw records and merge samples into one callset.

Validation is the single funnel for both input shapes. Merging rules:

* records with the same normalized :class:`VariantKey` are aggregated into one
  callset entry; carriers are de-duplicated, so a repeated/equivalent writing
  never inflates the sample count;
* within one sample, repeated identical calls are treated as supporting
  evidence, not conflicts;
* semantic conflicts are *marked*, never resolved by picking one record:
  HET vs HOM for the same variant (``GENOTYPE_CONFLICT``), or HOM_ALT vs
  homozygous-reference for the same anchored site (``REF_CONFLICT``);
* two different ALT alleles for the same sample/site are flagged via
  ``SAMPLE_CALL_CONFLICT`` on both entries (each ALT is real evidence, but the
  sample's call for the site is ambiguous);
* output ordering is deterministic.
"""
from __future__ import annotations

import re

from .models import (
    CallsetEntry,
    CallsetResult,
    Conflict,
    ConflictType,
    Genotype,
    Location,
    RawRecord,
    Rejection,
    RejectReason,
    VariantKey,
    VariantRecord,
)
from .normalization import (
    build_variant_key,
    canonicalize_genotype,
    is_valid_allele,
    normalize_chrom,
)

# genotype rank for choosing what a sample *calls* (strongest ALT wins)
_RANK = {Genotype.HOM: 4, Genotype.HET: 3, Genotype.MISSING: 2, Genotype.REF: 1}


def _reject(reason: RejectReason, location: Location, detail: str,
            sample: str | None = None) -> Rejection:
    return Rejection(reason=reason, location=location, detail=detail, sample=sample)


def validate_and_normalize(
    raws: list[RawRecord],
) -> tuple[list[VariantRecord], list[Rejection]]:
    """Semantically validate raw records; return accepted + rejected lists.

    Every raw record represents one concrete sample cell, so each semantic
    rejection is attributed to that sample; this lets an incremental update
    replace one sample's rejects along with its calls.
    """
    accepted: list[VariantRecord] = []
    rejected: list[Rejection] = []

    for raw in raws:
        loc = raw.location
        # chromosome
        try:
            chrom = normalize_chrom(raw.chrom)
        except ValueError:
            rejected.append(_reject(RejectReason.BAD_CHROM, loc, f"染色体非法: {raw.chrom!r}",
                                   raw.sample))
            continue
        # position
        pos_text = raw.pos_raw.strip()
        if not re.fullmatch(r"\d+", pos_text):
            rejected.append(_reject(RejectReason.BAD_POSITION, loc,
                                    f"POS 不是非负整数: {raw.pos_raw!r}", raw.sample))
            continue
        pos = int(pos_text)
        if pos < 1:
            rejected.append(_reject(RejectReason.BAD_POSITION, loc,
                                    f"POS 必须 >= 1: {pos}", raw.sample))
            continue
        # alleles
        ref, alt = raw.ref, raw.alt
        if not is_valid_allele(ref) or not is_valid_allele(alt):
            rejected.append(
                _reject(RejectReason.BAD_ALLELE, loc,
                        f"REF/ALT 含非法碱基或为空: {ref!r}/{alt!r}", raw.sample)
            )
            continue
        if ref == alt:
            rejected.append(
                _reject(RejectReason.BAD_ALLELE, loc,
                        f"REF 与 ALT 相同，不构成变异: {ref}", raw.sample)
            )
            continue
        # genotype
        try:
            genotype = canonicalize_genotype(raw.genotype)
        except ValueError:
            rejected.append(
                _reject(RejectReason.BAD_GENOTYPE, loc,
                        f"基因型无法识别: {raw.genotype!r}", raw.sample)
            )
            continue
        # normalization (left alignment); guard against any coordinate error
        try:
            key = build_variant_key(chrom, pos, ref, alt)
        except ValueError as exc:
            rejected.append(_reject(RejectReason.MALFORMED_LINE, loc, str(exc), raw.sample))
            continue
        accepted.append(VariantRecord(key=key, sample=raw.sample, genotype=genotype, location=loc))

    return accepted, rejected


def _site_key(key: VariantKey) -> tuple[str, int]:
    """Anchor site identity used for cross-ALT reference conflicts."""
    return (key.chrom, key.pos)


def _chrom_sort_token(chrom: str) -> tuple[int, int | str]:
    body = chrom[3:] if chrom.startswith("CHR") else chrom
    if body.isdigit():
        return (0, int(body))
    special = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if body.upper() in special:
        return (1, special[body.upper()])
    return (2, body)


def _sort_key(key: VariantKey):
    return (_chrom_sort_token(key.chrom), key.pos, key.ref, key.alt)


def _loc_key(loc: Location):
    return (loc.block, loc.line_no, loc.raw)


def _conflict_key(c: Conflict):
    return (
        c.conflict_type.value,
        c.variant.render(),
        c.sample or "",
        c.detail,
        tuple(_loc_key(l) for l in c.locations),
    )


def merge_callset(
    content_hash: str,
    accepted: list[VariantRecord],
    rejected: list[Rejection],
) -> CallsetResult:
    """Aggregate normalized records into a deterministic, conflict-annotated callset."""
    # key -> sample -> {rank, genotype, locations}
    groups: dict[VariantKey, dict[str, dict]] = {}
    # site -> sample -> set of ALT keys the sample has an ALT call for
    site_alts: dict[tuple[str, int], dict[str, set[VariantKey]]] = {}
    # site -> sample -> hom-ref locations
    site_ref: dict[tuple[str, int], dict[str, list[Location]]] = {}

    for rec in accepted:
        per_sample = groups.setdefault(rec.key, {}).setdefault(
            rec.sample, {"rank": 0, "genotype": Genotype.MISSING, "locs": [], "genotypes": set()}
        )
        per_sample["locs"].append(rec.location)
        per_sample["genotypes"].add(rec.genotype)
        if _RANK[rec.genotype] > per_sample["rank"]:
            per_sample["rank"] = _RANK[rec.genotype]
            per_sample["genotype"] = rec.genotype

        site = _site_key(rec.key)
        if rec.genotype in (Genotype.HET, Genotype.HOM):
            site_alts.setdefault(site, {}).setdefault(rec.sample, set()).add(rec.key)
        elif rec.genotype is Genotype.REF:
            site_ref.setdefault(site, {}).setdefault(rec.sample, []).append(rec.location)

    top_conflicts: dict[tuple, Conflict] = {}

    def add_conflict(c: Conflict) -> None:
        top_conflicts.setdefault(_conflict_key(c), c)

    entries: list[CallsetEntry] = []

    for key in sorted(groups, key=_sort_key):
        per_sample = groups[key]
        entry_conflicts: list[Conflict] = []
        genotype_by_sample: dict[str, Genotype] = {}
        carriers: list[str] = []
        support: list[Location] = []

        for sample in sorted(per_sample):
            info = per_sample[sample]
            genotype: Genotype = info["genotype"]
            locs: list[Location] = info["locs"]
            support.extend(locs)

            distinct: set[Genotype] = info["genotypes"]
            # GENOTYPE_CONFLICT: both HET and HOM observed for the same variant.
            if Genotype.HET in distinct and Genotype.HOM in distinct:
                c = Conflict(
                    conflict_type=ConflictType.GENOTYPE_CONFLICT,
                    variant=key,
                    sample=sample,
                    detail=(
                        f"样本 {sample} 在 {key.render()} 同时出现杂合(0/1)与纯合(1/1)，"
                        "不选择任何一条覆盖另一条"
                    ),
                    locations=tuple(sorted(locs, key=_loc_key)),
                )
                entry_conflicts.append(c)
                add_conflict(c)
                genotype = genotype  # keep strongest but flag the ambiguity
            # REF_CONFLICT: HOM_ALT here but a homozygous-reference call at same site.
            if genotype is Genotype.HOM:
                ref_locs = site_ref.get(_site_key(key), {}).get(sample, [])
                if ref_locs:
                    c = Conflict(
                        conflict_type=ConflictType.REF_CONFLICT,
                        variant=key,
                        sample=sample,
                        detail=(
                            f"样本 {sample} 在 {key.chrom}:{key.pos} 同时被记为纯合 ALT "
                            f"({key.render()}) 与纯合参考(0/0)"
                        ),
                        locations=tuple(sorted(locs + ref_locs, key=_loc_key)),
                    )
                    entry_conflicts.append(c)
                    add_conflict(c)

            genotype_by_sample[sample] = genotype
            if genotype in (Genotype.HET, Genotype.HOM):
                carriers.append(sample)

        # SAMPLE_CALL_CONFLICT: same sample has ALT calls for distinct ALTs at site.
        site = _site_key(key)
        for sample, alt_keys in site_alts.get(site, {}).items():
            if sample in per_sample and len(alt_keys) > 1:
                others = sorted(alt_keys - {key}, key=_sort_key)
                if others:
                    related_locs = [
                        r.location
                        for r in accepted
                        if _site_key(r.key) == site and r.sample == sample
                    ]
                    c = Conflict(
                        conflict_type=ConflictType.SAMPLE_CALL_CONFLICT,
                        variant=key,
                        sample=sample,
                        detail=(
                            f"样本 {sample} 在 {key.chrom}:{key.pos} 对不同 ALT 等位基因"
                            f"({key.alt} 与 {','.join(o.alt for o in others)}) 均有 ALT 调用"
                        ),
                        locations=tuple(sorted(related_locs, key=_loc_key)),
                    )
                    entry_conflicts.append(c)
                    add_conflict(c)

        entries.append(
            CallsetEntry(
                variant=key,
                samples=tuple(sorted(set(carriers))),
                genotype_by_sample=dict(sorted(genotype_by_sample.items())),
                support=tuple(sorted(support, key=_loc_key)),
                conflicts=tuple(sorted(entry_conflicts, key=_conflict_key)),
            )
        )

    ordered_conflicts = [top_conflicts[k] for k in sorted(top_conflicts)]
    return CallsetResult(
        content_hash=content_hash,
        entries=entries,
        rejected=sorted(rejected, key=lambda r: _loc_key(r.location)),
        conflicts=ordered_conflicts,
        accepted_count=len(accepted),
    )
