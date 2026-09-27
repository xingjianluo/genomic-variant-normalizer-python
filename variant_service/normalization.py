"""变异文本归一化。

支持两种输入形态：

1. ``single`` 单样本形态（每行一条，制表符分隔）::

       sample_id <TAB> chrom <TAB> pos <TAB> ref <TAB> alt [<TAB> genotype>]

2. ``vcf`` 多样本形态：标准 VCF 文本，``#CHROM ... FORMAT <samples...>`` 表头
   之后每行一个位点、每个样本列展开成一条 call；无样本列时整行作为一条
   未知样本（UNKNOWN）call。

归一规则（等价写法收敛）：
* 染色体大小写/前缀归一：``1``、``CHR1`` -> ``chr1``；性染色体/线粒体等价写法
  （``chrX``/``chrX``，``chrM`` -> ``chrMT``）；无法识别的染色体拒绝。
* 等位基因一律大写，只允许 ACGTN 字符，拒绝符号等位（如 ``<DEL>``）。
* 无参考基因组条件下做有界左移（trimming）：同时去掉 REF/ALT 公共后缀，
  再循环去掉公共前缀，保留 1 个锚定碱基，从而把 ``A|AA`` 与 ``CA|CAA`` 这类
  带锚定差异的等价插入写法收敛为同一表示。
* 基因型接受 ``0/1``、``0|1``、``1/1``、``./.`` 以及 ``HET``/``HOM_ALT`` 等
  规范名；多 ALT 位点的 GT 按等位索引展开为多条 call。
"""
from __future__ import annotations

from .models import Genotype, RejectReason, RejectedRecord, VariantCall
from .models import ParseResult

VALID_ALLELE = set("ACGTN")

# 展示与排序用的规范染色体顺序
CANONICAL_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrMT"]
_CHROM_ORDER = {c: i for i, c in enumerate(CANONICAL_CHROMS)}

_CHROM_ALIAS = {"M": "chrMT", "MT": "chrMT", "X": "chrX", "Y": "chrY"}

_DEFAULT_SAMPLE = "UNKNOWN"


def normalize_chrom(raw: str) -> str | None:
    """染色体名归一。无法识别返回 None。"""
    s = raw.strip().upper()
    if s.startswith("CHR"):
        s = s[3:]
    if s.isdigit() and 1 <= int(s) <= 22:
        return f"chr{s}"
    if s in _CHROM_ALIAS:
        return _CHROM_ALIAS[s]
    return None


def chrom_order(chrom: str) -> int:
    """规范染色体的排序权重；无法识别的染色体排在最后（按字典序稳定）。"""
    return _CHROM_ORDER.get(chrom, len(CANONICAL_CHROMS))


def normalize_allele(ref: str, alt: str, pos: int) -> tuple[int, str, str] | None:
    """等位基因校验 + 有界左移。

    返回 ``(pos, ref, alt)`` 归一结果；任一碱基非法返回 None。
    """
    ref = ref.upper()
    alt = alt.upper()
    if not ref or not alt or ref == alt:
        return None
    if not (set(ref) <= VALID_ALLELE and set(alt) <= VALID_ALLELE):
        return None

    # 1) 去掉公共前缀（保留至少一个锚定碱基），左移 POS
    trim = 0
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        trim += 1
    # 2) 去掉公共后缀（保留至少一个锚定碱基）
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos + trim, ref, alt


_GT_NAME = {
    "HOM_REF": Genotype.HOM_REF,
    "HET": Genotype.HET,
    "HOM_ALT": Genotype.HOM_ALT,
    "UNKNOWN": Genotype.UNKNOWN,
}


def parse_genotype(raw: str | None) -> Genotype | None:
    """解析基因型文本。返回 None 表示无基因型（按 UNKNOWN 处理）。

    非法数字基因型抛出 :class:`ValueError`，调用方据此分类 INVALID_GENOTYPE。
    """
    if raw is None or raw.strip() == "":
        return Genotype.UNKNOWN
    s = raw.strip().upper()
    if s in _GT_NAME:
        return _GT_NAME[s]
    if s in {".", "./.", ".|.", "././.", ".|.|."}:
        return Genotype.UNKNOWN

    # 数字 GT，如 0/1、1|0、1/1
    if any(sep in s for sep in ("/", "|")):
        tokens = s.replace("|", "/").split("/")
        if any(t in ("", ".") for t in tokens):
            # 部分缺失（如 ./1）：作为 UNKNOWN 处理，不强行判型
            return Genotype.UNKNOWN
        try:
            idxs = [int(t) for t in tokens]
        except ValueError:
            raise ValueError(f"非法基因型: {raw!r}")
        if any(i < 0 for i in idxs):
            raise ValueError(f"非法基因型: {raw!r}")
        non_ref = [i for i in idxs if i != 0]
        if not non_ref:
            return Genotype.HOM_REF
        if len(set(non_ref)) == 1 and len(idxs) > 1 and len(non_ref) == len(idxs):
            return Genotype.HOM_ALT
        return Genotype.HET
    raise ValueError(f"无法识别的基因型: {raw!r}")


def gt_for_indices(idxs: list[int]) -> Genotype:
    """根据 VCF GT 等位索引（均已解析为非负整数）判定基因型。"""
    non_ref = [i for i in idxs if i != 0]
    if not non_ref:
        return Genotype.HOM_REF
    if len(set(non_ref)) == 1 and len(idxs) > 1 and len(non_ref) == len(idxs):
        return Genotype.HOM_ALT
    return Genotype.HET


def _reject(source: str, line: int, reason: RejectReason, detail: str, raw: str) -> RejectedRecord:
    return RejectedRecord(source=source, line=line, reason=reason, detail=detail, raw=raw)


def _build_call(
    result: ParseResult,
    source: str,
    line_no: int,
    sample_id: str,
    chrom_raw: str,
    pos_raw: str,
    ref_raw: str,
    alt_raw: str,
    gt_raw: str | None,
    raw_line: str,
    alt_index: int | None = None,
) -> None:
    """校验并归一一条记录，写入 result（成功为 call，失败为 rejection）。"""
    if not chrom_raw.strip() or not pos_raw.strip() or not ref_raw.strip() or not alt_raw.strip():
        result.rejected.append(
            _reject(source, line_no, RejectReason.MISSING_FIELD, "chrom/pos/ref/alt 不允许为空", raw_line)
        )
        return
    chrom = normalize_chrom(chrom_raw)
    if chrom is None:
        result.rejected.append(
            _reject(source, line_no, RejectReason.INVALID_CHROM, f"无法识别的染色体: {chrom_raw!r}", raw_line)
        )
        return
    try:
        pos = int(pos_raw.strip())
    except ValueError:
        result.rejected.append(
            _reject(source, line_no, RejectReason.INVALID_POSITION, f"位置不是整数: {pos_raw!r}", raw_line)
        )
        return
    if pos < 1:
        result.rejected.append(
            _reject(source, line_no, RejectReason.INVALID_POSITION, f"位置必须 >=1: {pos}", raw_line)
        )
        return
    norm = normalize_allele(ref_raw.strip(), alt_raw.strip(), pos)
    if norm is None:
        result.rejected.append(
            _reject(source, line_no, RejectReason.INVALID_ALLELE,
                    f"非法等位基因 ref={ref_raw!r} alt={alt_raw!r}（仅支持 ACGTN 且需不同）", raw_line)
        )
        return
    npos, nref, nalt = norm

    # 基因型
    if gt_raw is None:
        genotype = Genotype.UNKNOWN
    else:
        try:
            genotype = parse_genotype(gt_raw)
        except ValueError as exc:
            result.rejected.append(
                _reject(source, line_no, RejectReason.INVALID_GENOTYPE, str(exc), raw_line)
            )
            return
        # 数字 GT 的等位索引必须在 ALT 数量范围内
        if alt_index is not None and genotype is not Genotype.UNKNOWN:
            s = gt_raw.strip().upper().replace("|", "/")
            if "/" in s:
                tokens = [t for t in s.split("/") if t not in ("", ".")]
                try:
                    idxs = [int(t) for t in tokens]
                    if any(i > alt_index for i in idxs if i != 0):
                        result.rejected.append(
                            _reject(source, line_no, RejectReason.INVALID_GENOTYPE,
                                    f"GT 等位索引超出 ALT 数量: {gt_raw!r}", raw_line)
                        )
                        return
                except ValueError:
                    pass  # 已在 parse_genotype 中拒绝

    result.calls.append(
        VariantCall(
            sample_id=sample_id,
            chrom=chrom,
            pos=npos,
            ref=nref,
            alt=nalt,
            genotype=genotype,
            source_line=line_no,
        )
    )


def parse_single(sample_id: str, text: str) -> ParseResult:
    """解析单样本形态文本。``sample_id`` 为空时用载荷首列作为样本名。"""
    result = ParseResult()
    source = sample_id or "single-input"
    payload_has_sample = sample_id == ""

    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if payload_has_sample:
            if len(parts) not in (5, 6):
                result.rejected.append(
                    _reject(source, line_no, RejectReason.WRONG_COLUMN_COUNT,
                            f"期望 5 或 6 列（含样本列），实际 {len(parts)} 列", raw_line)
                )
                continue
            sid, chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            gt = parts[5] if len(parts) == 6 else None
            if not sid.strip():
                result.rejected.append(
                    _reject(source, line_no, RejectReason.MISSING_FIELD, "样本名为空", raw_line)
                )
                continue
        else:
            if len(parts) not in (4, 5):
                result.rejected.append(
                    _reject(source, line_no, RejectReason.WRONG_COLUMN_COUNT,
                            f"期望 4 或 5 列（chrom/pos/ref/alt[+gt]），实际 {len(parts)} 列", raw_line)
                )
                continue
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            gt = parts[4] if len(parts) == 5 else None
            sid = sample_id
        _build_call(result, sid or source, line_no, sid, chrom, pos, ref, alt, gt, raw_line)
    return result


def parse_vcf(text: str) -> ParseResult:
    """解析多样本 VCF 文本。"""
    result = ParseResult()
    source = "vcf"
    sample_names: list[str] | None = None

    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        if raw_line.startswith("##"):
            continue
        if raw_line.startswith("#"):
            header = raw_line.lstrip("#").split("\t")
            if len(header) < 8 or header[0].strip().upper() != "CHROM":
                result.rejected.append(
                    _reject(source, line_no, RejectReason.MALFORMED_LINE,
                            "VCF 表头必须以 #CHROM 开头且至少含 8 个固定列", raw_line)
                )
                continue
            names = [h.strip() for h in header[9:]]
            if len(names) != len(set(names)):
                dupes = sorted({n for n in names if names.count(n) > 1})
                result.rejected.append(
                    _reject(source, line_no, RejectReason.DUPLICATE_SAMPLE_HEADER,
                            f"样本列重名: {dupes}", raw_line)
                )
                continue
            sample_names = names
            continue

        line = raw_line.strip()
        if not line:
            continue
        if sample_names is None:
            result.rejected.append(
                _reject(source, line_no, RejectReason.MALFORMED_LINE,
                        "数据行之前缺少 #CHROM 表头", raw_line)
            )
            continue

        parts = line.split("\t")
        if len(parts) < 8:
            result.rejected.append(
                _reject(source, line_no, RejectReason.WRONG_COLUMN_COUNT,
                        f"VCF 数据行至少需要 8 列，实际 {len(parts)} 列", raw_line)
            )
            continue
        chrom, pos, _vid, ref, alt_field, _qual, _filt, info = parts[:8]
        alts = [a.strip() for a in alt_field.split(",")]

        if len(parts) == 8:
            # 无样本列：整行作为一条 UNKNOWN 样本 / UNKNOWN 基因型 call（每个 ALT）
            for alt in alts:
                _build_call(result, source, line_no, _DEFAULT_SAMPLE,
                            chrom, pos, ref, alt, None, raw_line)
            continue

        if len(parts) == 9:
            result.rejected.append(
                _reject(source, line_no, RejectReason.WRONG_COLUMN_COUNT,
                        "存在 FORMAT 列但缺少样本值列（第 9 列后无样本）", raw_line)
            )
            continue

        format_fields = parts[8].split(":")
        sample_cells = parts[9:]
        if len(sample_cells) != len(sample_names):
            result.rejected.append(
                _reject(source, line_no, RejectReason.SAMPLE_COLUMN_MISMATCH,
                        f"样本值列数 {len(sample_cells)} 与表头样本数 {len(sample_names)} 不一致",
                        raw_line)
            )
            continue

        try:
            gt_pos = format_fields.index("GT")
        except ValueError:
            gt_pos = -1

        for sid, cell in zip(sample_names, sample_cells):
            vals = cell.split(":")
            gt_raw = vals[gt_pos] if gt_pos >= 0 and gt_pos < len(vals) else None
            gt_text = gt_raw.strip() if gt_raw is not None else ""
            if gt_text in ("", ".", "./.", ".|."):
                genotype = Genotype.UNKNOWN
                numeric_idxs: list[int] | None = None
            else:
                tokens = gt_text.replace("|", "/").split("/")
                if all(t.isdigit() for t in tokens):
                    numeric_idxs = [int(t) for t in tokens]
                    if any(i > len(alts) for i in numeric_idxs if i != 0):
                        result.rejected.append(
                            _reject(sid, line_no, RejectReason.INVALID_GENOTYPE,
                                    f"GT 等位索引超出 ALT 数量: {gt_text!r}", raw_line)
                        )
                        continue
                    genotype = gt_for_indices(numeric_idxs)
                else:
                    try:
                        genotype = parse_genotype(gt_raw)
                    except ValueError as exc:
                        result.rejected.append(
                            _reject(sid, line_no, RejectReason.INVALID_GENOTYPE, str(exc), raw_line)
                        )
                        continue
                    numeric_idxs = None

            # 该样本在此位点涉及哪些 ALT：
            # - HET/HOM_ALT：按 GT 中的非零索引展开
            # - HOM_REF / UNKNOWN：对每个 ALT 各发一条证据 call
            if genotype in (Genotype.HOM_REF, Genotype.UNKNOWN):
                alt_indices = list(range(len(alts)))
            else:
                alt_indices = sorted({i - 1 for i in numeric_idxs if i >= 1})

            for ai in alt_indices:
                _build_call(result, sid, line_no, sid,
                            chrom, pos, ref, alts[ai],
                            genotype.value, raw_line, alt_index=len(alts))
    return result


def parse_input(sample_id: str, fmt: str, text: str) -> ParseResult:
    """根据格式分发解析；空文本给 EMPTY_INPUT 拒绝。"""
    if text is None or not text.strip():
        reason_source = sample_id or fmt
        result = ParseResult()
        result.rejected.append(
            RejectedRecord(source=reason_source, line=0, reason=RejectReason.EMPTY_INPUT,
                           detail="输入为空或仅含空白", raw="")
        )
        return result
    if fmt == "vcf":
        return parse_vcf(text)
    return parse_single(sample_id, text)
