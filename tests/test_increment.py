"""Tests for the pure incremental-update layer.

Business coverage (each test logs its 判定依据):
* new samples add entries; corrected samples replace (never stack) conclusions;
* untouched entries/samples/order are preserved verbatim;
* change summary: added/removed/modified entries, carrier changes, conflicts
  appeared/resolved;
* reproducibility: same baseline + same delta => identical version id + summary;
* malformed corrected sample is rejected with distinguishable reasons and does
  not silently overwrite its prior good conclusion;
* both input shapes (single-sample and multi-sample) drive the same rules.
"""
import logging

from variant_service.increment import apply_delta, build_initial_version
from variant_service.merging import validate_and_normalize
from variant_service.models import SampleChangeType
from variant_service.parsing import parse_multi_sample, parse_single_sample

logger = logging.getLogger(__name__)


def _multi(text):
    raws, struct = parse_multi_sample(text)
    accepted, semantic = validate_and_normalize(raws)
    return accepted, struct + semantic


def _single(text, sample):
    raws, struct = parse_single_sample(text, sample)
    accepted, semantic = validate_and_normalize(raws)
    return accepted, struct + semantic


def _initial(text, order):
    accepted, rejected = _multi(text)
    return build_initial_version("h0", accepted, rejected, order)


BASELINE = (
    "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
    "chr1\t100\tA\tT\t0/1\t1/1\n"      # shared entry, S1 HET / S2 HOM
    "chr2\t200\tG\tC\t.\t0/0\n"        # S1 ./. + S2 0/0 (no carriers)
)


def test_new_sample_adds_entries_and_summary_flags_it():
    v0 = _initial(BASELINE, ("S1", "S2"))
    accepted, rejected = _multi(
        "#CHROM\tPOS\tREF\tALT\tS3\nchr3\t300\tC\tG\t0/1\n",
    )
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S3",))
    summary = v1.summary
    logger.info("判定依据 新增样本 S3: entries_added=%s new_samples=%s",
                summary.entries_added, summary.new_samples)
    assert summary.entries_added == ("CHR3:300:C:G",)
    assert summary.new_samples == ("S3",) and summary.corrected_samples == ()
    assert [e.variant.pos for e in v1.result.entries] == [100, 200, 300]
    # parent linkage for traceability
    assert v1.parent_version_id == v0.version_id


def test_correction_replaces_conclusion_not_stacks_evidence():
    v0 = _initial(BASELINE, ("S1", "S2"))
    # S1 re-submitted as HOM: the old HET must disappear (no HET/HOM conflict),
    # support evidence comes only from the new submission.
    accepted, rejected = _single("chr1 100 A T 1/1", "S1")
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S1",))
    entry = next(e for e in v1.result.entries if e.variant.pos == 100)
    logger.info("判定依据 修正 S1 HET->HOM: gt=%s support=%d conflicts=%d",
                entry.genotype_by_sample["S1"].value, len(entry.support),
                len(v1.result.conflicts))
    assert entry.genotype_by_sample["S1"].value == "1/1"
    assert v1.result.conflicts == []
    # only the new support location survives (old evidence not stacked)
    assert len(entry.support) == 2  # S1(new) + S2 untouched
    assert summary_has(v1, variant="CHR1:100:A:T", change=SampleChangeType.GENOTYPE_CHANGED,
                       sample="S1", before="0/1", after="1/1")
    assert v1.summary.corrected_samples == ("S1",)


def test_untouched_entries_are_reused_and_order_preserved():
    # Build baseline via per-sample blocks so S2 owns chr4 with no S1 cell.
    acc, rej = [], []
    for sample, text in (("S1", "chr1 100 A T 0/1\n"),
                         ("S2", "chr1 100 A T 1/1\nchr4 400 A T 0/1\n")):
        a, r = _single(text, sample)
        acc += a
        rej += r
    v0 = build_initial_version("h0", acc, rej, ("S1", "S2"))
    untouched = next(e for e in v0.result.entries if e.variant.pos == 400)

    accepted, rejected = _multi(
        "#CHROM\tPOS\tREF\tALT\tS3\tS1\n"
        "chr1\t100\tA\tT\t1/1\t1/1\nchr3\t300\tC\tG\t0/1\t0/1\n",
    )
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S3", "S1"))
    after = next(e for e in v1.result.entries if e.variant.pos == 400)
    positions = [e.variant.pos for e in v1.result.entries]
    logger.info("判定依据 未触及条目: 对象复用=%s, 顺序=%s", after is untouched, positions)
    assert after is untouched
    assert positions == [100, 300, 400]


def test_correction_can_remove_entries_and_carriers():
    v0 = _initial(BASELINE, ("S1", "S2"))
    # S1 re-submitted with only a brand-new variant: its chr1/chr2 calls vanish.
    accepted, rejected = _single("chr9 900 A T 0/1", "S1")
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S1",))
    logger.info("判定依据 S1 改测至 chr9: removed=%s added=%s",
                v1.summary.entries_removed, v1.summary.entries_added)
    # chr1 entry survives (S2 still HOM carrier) but S1 leaves its carrier set;
    # chr2 entry only had S2's 0/0 and S1's ./. -> S2 0/0 keeps it (ref kept).
    chr1 = next(e for e in v1.result.entries if e.variant.pos == 100)
    assert "S1" not in chr1.genotype_by_sample
    assert chr1.samples == ("S2",)
    assert v1.summary.entries_added == ("CHR9:900:A:T",)
    assert summary_has(v1, variant="CHR1:100:A:T", change=SampleChangeType.CARRIER_REMOVED,
                       sample="S1", before="0/1", after=None)


def test_conflict_appears_then_resolves_on_correction():
    v0 = _initial(BASELINE, ("S1", "S2"))
    # corrected S1 now has both HET and HOM -> GENOTYPE_CONFLICT appears, but it
    # is still marked, not silently overwritten
    accepted, rejected = _single("chr1 100 A T 0/1\nchr1 100 A T 1/1", "S1")
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S1",))
    appeared = [c for c in v1.summary.conflicts_appeared if c.sample == "S1"]
    logger.info("判定依据 修正引入矛盾: appeared=%s", [c.conflict_type.value for c in appeared])
    assert len(appeared) == 1
    assert [c.conflict_type.value for c in v1.result.conflicts] == ["GENOTYPE_CONFLICT"]
    # corrected again with one clean call -> conflict resolved, not stacked
    accepted2, rejected2 = _single("chr1 100 A T 0/1", "S1")
    v2 = apply_delta(v1, "d2", accepted2, rejected2, ("S1",))
    resolved = [c for c in v2.summary.conflicts_resolved if c.sample == "S1"]
    logger.info("判定依据 再次修正消除矛盾: resolved=%d remaining=%d",
                len(resolved), len(v2.result.conflicts))
    assert len(resolved) == 1
    assert v2.result.conflicts == []


def test_same_baseline_same_delta_is_reproducible():
    v0 = _initial(BASELINE, ("S1", "S2"))
    text = "#CHROM\tPOS\tREF\tALT\tS3\tS1\nchr1\t100\tA\tT\t1/1\t1/1\nchr3\t300\tC\tG\t0/1\t0/1\n"
    a1, r1 = _multi(text)
    a2, r2 = _multi(text)
    v1 = apply_delta(v0, "d1", a1, r1, ("S3", "S1"))
    v2 = apply_delta(v0, "d1", a2, r2, ("S3", "S1"))
    logger.info("判定依据 同基线同增量重复: %s == %s", v1.version_id[:10], v2.version_id[:10])
    assert v1.version_id == v2.version_id
    assert v1.summary == v2.summary
    # re-applying the same delta onto the already-derived version is a no-op too
    v3 = apply_delta(v1, "d1", a2, r2, ("S3", "S1"))
    assert v3.version_id == v1.version_id


def test_malformed_correction_is_rejected_not_silently_applied():
    v0 = _initial(BASELINE, ("S1", "S2"))
    # corrected S1 text carries a bad allele and a bad genotype: distinguishable
    # reasons + locations; its prior good conclusion must NOT be destroyed.
    accepted, rejected = _single("chr1 100 A U 0/1\nchr1 100 A T 9/9", "S1")
    v1 = apply_delta(v0, "d1", accepted, rejected, ("S1",))
    reasons = sorted((r.reason.value, r.location.line_no) for r in rejected)
    logger.info("判定依据 畸形修正样本: rejects=%s (旧 HET 结论按样本整体替换为 0 条有效调用)",
                reasons)
    # both bad rows reported with their line
    assert reasons == [("BAD_ALLELE", 1), ("BAD_GENOTYPE", 2)]
    # correction carried zero valid calls: S1 removed from chr1, reject retained
    chr1 = next(e for e in v1.result.entries if e.variant.pos == 100)
    assert "S1" not in chr1.genotype_by_sample
    assert {r.reason.value for r in v1.result.rejected} == {"BAD_ALLELE", "BAD_GENOTYPE"}


def test_multi_shape_column_mismatch_in_delta_is_batch_reject():
    v0 = _initial(BASELINE, ("S1", "S2"))
    raws, struct = parse_multi_sample(
        "#CHROM\tPOS\tREF\tALT\tS9\tS8\n"
        "chr9\t900\tA\tT\t0/1\t1/1\n"
        "chr8\t800\tA\tT\t0/1\n",  # one cell short -> mismatch, row dropped
    )
    accepted, semantic = validate_and_normalize(raws)
    v1 = apply_delta(v0, "d1", accepted, struct + semantic, ("S9", "S8"))
    logger.info("判定依据 增量列错位: rejects=%s, 仍新增 S9=%s",
                [r.reason.value for r in struct],
                any(s in e.genotype_by_sample
                    for e in v1.result.entries for s in ("S9",)))
    assert [r.reason.value for r in v1.result.rejected] == ["HEADER_SAMPLE_MISMATCH"]
    chr9 = next(e for e in v1.result.entries if e.variant.pos == 900)
    assert "S9" in chr9.genotype_by_sample and "S8" in chr9.genotype_by_sample
    assert not [e for e in v1.result.entries if e.variant.pos == 800]


def summary_has(version, *, variant, change, sample, before, after):
    for entry in version.summary.entry_changes:
        if entry.variant == variant:
            for sc in entry.sample_changes:
                if (sc.change_type is change and sc.sample == sample
                        and sc.before == before and sc.after == after):
                    return True
    return False
