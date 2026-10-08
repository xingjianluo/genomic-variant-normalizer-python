"""Parse single-sample and multi-sample variant text into the shared model.

Parsing is *structural* only: it splits rows, locates columns and preserves raw
field text. Semantic validation (position integer, allele alphabet, genotype
vocabulary, REF!=ALT) happens later in :mod:`variant_service.merging`, so the
two input shapes converge on exactly one validation path.

Supported shapes
----------------
Single sample (one block, repeated across requests for different samples)::

    # comment
    chr1  100  A    T     0/1
    chr2  200  ACG  AT

Multi-sample, compact header::

    #CHROM  POS  REF  ALT  S1    S2
    chr1    100  A    T    0/1   1/1

Multi-sample, VCF-style header (fixed columns + FORMAT)::

    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2
    chr1   100 .  A   T   .    .      .    GT     0/1 1/1

A multi-sample block without a header is accepted: the sample count is inferred
from the first data row (samples named ``sample_1`` ...) and every later row
with a different column count is rejected as HEADER_SAMPLE_MISMATCH.
"""
from __future__ import annotations

from .models import Location, RawRecord, RejectReason, Rejection

_VCF_FIXED = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


def _loc(block: str, line_no: int, raw: str) -> Location:
    return Location(block=block, line_no=line_no, raw=raw)


def _reject(reason: RejectReason, location: Location, detail: str) -> Rejection:
    return Rejection(reason=reason, location=location, detail=detail)


def _split_row(line: str) -> list[str]:
    """Split a data row on tabs, falling back to whitespace; strip cells."""
    if "\t" in line:
        parts = line.split("\t")
    else:
        parts = line.split()
    return [p.strip() for p in parts]


def _raw_record(parts: list[str], gt: str, sample: str, location: Location) -> RawRecord:
    return RawRecord(
        chrom=parts[0],
        pos_raw=parts[1],
        ref=parts[2].upper(),
        alt=parts[3].upper(),
        genotype=gt.strip(),
        sample=sample,
        location=location,
    )


def parse_single_sample(text: str, sample: str) -> tuple[list[RawRecord], list[Rejection]]:
    """Parse one single-sample block into raw records and structural rejections."""
    block = f"sample:{sample}"
    raws: list[RawRecord] = []
    rejects: list[Rejection] = []

    for line_no, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        loc = _loc(block, line_no, line.rstrip("\r"))
        parts = _split_row(stripped)
        if len(parts) < 4:
            rejects.append(
                _reject(RejectReason.MISSING_FIELDS, loc,
                        f"需要至少 CHROM/POS/REF/ALT 四列，实际 {len(parts)} 列")
            )
            continue
        if len(parts) > 5:
            rejects.append(
                _reject(RejectReason.MALFORMED_LINE, loc,
                        f"单样本行至多 5 列(含可选 GT)，实际 {len(parts)} 列")
            )
            continue
        gt = parts[4] if len(parts) == 5 else ""
        try:
            raws.append(_raw_record(parts, gt, sample, loc))
        except IndexError:
            rejects.append(_reject(RejectReason.MALFORMED_LINE, loc, "字段索引失败"))

    return raws, rejects


def _dedupe_names(names: list[str], header_loc: Location) -> tuple[list[str], list[Rejection]]:
    seen: dict[str, int] = {}
    out: list[str] = []
    rejects: list[Rejection] = []
    flagged: set[str] = set()
    for name in names:
        if name in seen:
            seen[name] += 1
            unique = f"{name}_{seen[name]}"
            out.append(unique)
            if name not in flagged:
                flagged.add(name)
                rejects.append(
                    _reject(RejectReason.DUPLICATE_SAMPLE, header_loc,
                            f"表头中样本名重复: {name!r}，后续列改记为 {unique!r}")
                )
        else:
            seen[name] = 1
            out.append(name)
    return out, rejects


def _parse_header(line: str) -> tuple[list[str], bool, int]:
    """Return (sample_names, is_vcf_style, gt_index_in_format).

    Compact header has 4 fixed columns then sample names. VCF header has the 8
    standard fixed columns, an optional FORMAT column, then sample names.
    """
    cols = _split_row(line)
    upper = [c.upper() for c in cols]
    if len(cols) >= 8 and tuple(upper[:8]) == _VCF_FIXED:
        gt_idx = 0
        if len(cols) > 8 and upper[8] == "FORMAT":
            format_fields = [f.strip().upper() for f in cols[8].split(":")]
            gt_idx = format_fields.index("GT") if "GT" in format_fields else -1
            names = cols[9:]
        else:
            names = cols[8:]
        return names, True, gt_idx
    # compact: #CHROM POS REF ALT samples...
    return cols[4:], False, 0


def parse_multi_sample(text: str) -> tuple[list[RawRecord], list[Rejection]]:
    """Parse a multi-sample block (optional header) into raw records/rejections."""
    block = "multisample"
    raws: list[RawRecord] = []
    rejects: list[Rejection] = []

    names: list[str] | None = None
    vcf_style = False
    gt_idx = 0
    inferred_columns: int | None = None  # total columns incl. 4 fixed

    for line_no, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("##"):
            continue
        loc = _loc(block, line_no, line.rstrip("\r"))

        if stripped.startswith("#"):
            header_names, vcf_style, gt_idx = _parse_header(stripped[1:])
            if not header_names:
                # A header without any sample column is not a multi-sample block.
                rejects.append(
                    _reject(RejectReason.MALFORMED_LINE, loc, "表头缺少样本列")
                )
                continue
            names, dup_rejects = _dedupe_names(header_names, loc)
            rejects.extend(dup_rejects)
            continue

        parts = _split_row(stripped)
        fixed = 9 if vcf_style else 4
        if len(parts) < fixed:
            rejects.append(
                _reject(RejectReason.MISSING_FIELDS, loc,
                        f"需要至少 {fixed} 个固定列，实际 {len(parts)} 列")
            )
            continue

        if names is None:
            # No header seen: infer sample count from this first data row.
            sample_count = len(parts) - 4
            if sample_count <= 0:
                rejects.append(
                    _reject(RejectReason.MISSING_FIELDS, loc, "无表头且数据行缺少样本列")
                )
                continue
            names = [f"sample_{i + 1}" for i in range(sample_count)]
            inferred_columns = len(parts)

        if vcf_style:
            # VCF style always has 8 fixed columns + FORMAT + one cell/sample.
            expected = 9 + len(names)
            if len(parts) != expected:
                rejects.append(
                    _reject(RejectReason.HEADER_SAMPLE_MISMATCH, loc,
                            f"样本列数不符: 表头 {len(names)} 个样本，期望 {expected} 列，"
                            f"实际 {len(parts)} 列")
                )
                continue
            format_fields = [f.strip().upper() for f in parts[8].split(":")] if parts[8] else []
            if "GT" not in format_fields:
                rejects.append(
                    _reject(RejectReason.MALFORMED_LINE, loc,
                            f"FORMAT={parts[8]!r} 缺少 GT 字段")
                )
                continue
            gt_idx = format_fields.index("GT")
            sample_cells = parts[9:]
            for name, cell in zip(names, sample_cells):
                gt_parts = cell.split(":")
                gt = gt_parts[gt_idx] if gt_idx < len(gt_parts) else ""
                raws.append(_raw_record(parts, gt, name, loc))
        else:
            expected = 4 + len(names)
            if inferred_columns is not None:
                expected = inferred_columns
            if len(parts) != expected:
                rejects.append(
                    _reject(RejectReason.HEADER_SAMPLE_MISMATCH, loc,
                            f"样本列数不符: 期望 {expected} 列（{len(names)} 个样本），"
                            f"实际 {len(parts)} 列")
                )
                continue
            sample_cells = parts[4:]
            for name, cell in zip(names, sample_cells):
                raws.append(_raw_record(parts, cell, name, loc))

    return raws, rejects


def declared_samples(text: str) -> tuple[str, ...]:
    """Return the ordered sample names a submission declares.

    Used by incremental updates to know which samples one batch carries (and
    therefore replaces). Mirrors :func:`parse_multi_sample` naming exactly:
    an explicit header (compact or VCF, with duplicate-name de-duplication),
    otherwise names inferred from the first data row's sample columns.
    Returns an empty tuple when the text has no usable header/data row.
    """
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("##"):
            continue
        if stripped.startswith("#"):
            header_names, _, _ = _parse_header(stripped[1:])
            if header_names:
                deduped, _ = _dedupe_names(header_names, _loc("multisample", 0, stripped))
                return tuple(deduped)
            continue
        parts = _split_row(stripped)
        if len(parts) > 4:
            return tuple(f"sample_{i + 1}" for i in range(len(parts) - 4))
        return ()
    return ()
