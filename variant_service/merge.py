"""多样本合并与冲突检测。

合并键为归一后的 ``(chrom, pos, ref, alt)`` 四元组；同键记录聚合成一个
callset 变异。规则：

* 去重：同一 ``(样本, 变异, 基因型)`` 的重复记录（含等价写法收敛后的重复）
  只计一次，不会重复计数。
* 基因型冲突（``GENOTYPE_CONFLICT``）：同一样本在同一变异上出现互相矛盾的
  具体基因型（HOM_REF 与带 ALT 的型别冲突；HET 与 HOM_ALT 冲突）。
  UNKNOWN 不与任何型别冲突。所有矛盾证据都保留在 genotype_summary 中，
  绝不任选一条覆盖另一条。
* 位点冲突（``SITE_CONFLICT``）：同一 ``(chrom, pos)`` 出现不同 REF 的记录。
  REF 相同 ALT 不同属于正常的多 ALT 位点，不算冲突。
* 输出按 (染色体规范顺序, pos, ref, alt) 稳定排序；同样输入必得同样输出。
"""
from __future__ import annotations

from collections import defaultdict

from .models import (
    Callset,
    CallsetVariant,
    Genotype,
    RejectedRecord,
    VariantCall,
)
from .normalization import chrom_order

# 定义具体基因型之间的互斥关系
_CONCRETE = {Genotype.HOM_REF, Genotype.HET, Genotype.HOM_ALT}


def variant_key(call: VariantCall) -> tuple[str, int, str, str]:
    return (call.chrom, call.pos, call.ref, call.alt)


def _genotypes_conflict(g1: Genotype, g2: Genotype) -> bool:
    if g1 is g2 or Genotype.UNKNOWN in (g1, g2):
        return False
    if g1 not in _CONCRETE or g2 not in _CONCRETE:
        return False
    return g1 is not g2


def merge_calls(
    callset_id: str,
    total_inputs: int,
    calls: list[VariantCall],
    rejected: list[RejectedRecord],
) -> Callset:
    # variant key -> sample -> set(genotypes)；用 set 天然去重
    grouped: dict[tuple[str, int, str, str], dict[str, set[Genotype]]] = defaultdict(
        lambda: defaultdict(set)
    )
    # 同一样本的同一 (sample,key,genotype) 重复只计一次
    seen: set[tuple[str, tuple, Genotype]] = set()
    accepted_calls = 0

    for call in calls:
        key = variant_key(call)
        sample_calls = grouped[key][call.sample_id]
        sig = (call.sample_id, key, call.genotype)
        if sig not in seen:
            seen.add(sig)
            accepted_calls += 1
        sample_calls.add(call.genotype)

    # 位点冲突检测：(chrom, pos) -> set(ref)
    site_refs: dict[tuple[str, int], set[str]] = defaultdict(set)
    for chrom, pos, ref, _alt in grouped:
        site_refs[(chrom, pos)].add(ref)
    conflicting_sites = {site for site, refs in site_refs.items() if len(refs) > 1}

    variants: list[CallsetVariant] = []
    for key, per_sample in grouped.items():
        chrom, pos, ref, alt = key
        genotype_summary: dict[str, list[str]] = {}
        samples: list[str] = []
        conflict_types: list[str] = []
        conflicts: list[dict] = []

        geno_conflict_samples: list[str] = []
        for sample_id in sorted(per_sample):
            gtypes = per_sample[sample_id]
            genotype_summary[sample_id] = sorted(g.value for g in gtypes)
            if any(g is not Genotype.UNKNOWN for g in gtypes):
                samples.append(sample_id)
            concrete = [g for g in gtypes if g in _CONCRETE]
            if len(set(concrete)) > 1:
                geno_conflict_samples.append(sample_id)
                conflicts.append({
                    "type": "GENOTYPE_CONFLICT",
                    "sample": sample_id,
                    "genotypes": sorted(g.value for g in gtypes),
                })

        if geno_conflict_samples:
            conflict_types.append("GENOTYPE_CONFLICT")

        if (chrom, pos) in conflicting_sites:
            conflict_types.append("SITE_CONFLICT")
            conflicts.append({
                "type": "SITE_CONFLICT",
                "chrom": chrom,
                "pos": pos,
                "refs": sorted(site_refs[(chrom, pos)]),
            })

        # 该变异的去重 call 数（样本×基因型组合数）
        call_count = sum(len(gts) for gts in per_sample.values())

        variants.append(
            CallsetVariant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                genotype_summary=dict(sorted(genotype_summary.items())),
                samples=samples,
                call_count=call_count,
                conflicting=bool(conflict_types),
                conflict_types=conflict_types,
                conflicts=conflicts,
            )
        )

    variants.sort(key=lambda v: (chrom_order(v.chrom), v.pos, v.ref, v.alt))

    rejected_sorted = sorted(rejected, key=lambda r: (r.source, r.line, r.reason.value))

    return Callset(
        callset_id=callset_id,
        total_inputs=total_inputs,
        accepted_calls=accepted_calls,
        rejected_records=rejected_sorted,
        variants=variants,
    )
