# 离线基因组变异记录处理服务

在原有最小 FastAPI 服务上扩展而成。把**多份样本**的变异文本归一化、合并成一份
可解释的 callset。全过程只使用本地数据与仓库自带能力（Python 标准库 + FastAPI），
**不依赖任何外部服务**（无网络注解、无 lift-over、无外部数据库）。

## 1. 支持的输入形态

两种形态最终统一到同一套记录模型（`variant_service/models.py`）：

### 1.1 单样本（`POST /normalize`，多次提交即可表达多样本）

每行 4~5 列，空白或制表符分隔；`#` 开头为注释，空行跳过：

```
# CHROM  POS  REF  ALT  [GT]
chr1     100  A    T    0/1
1        200  ACG  AT   hom
```

- 第 5 列基因型可省略（记为 `./.`，即“未给出”）。
- 请求需携带样本名 `sample`。

### 1.2 多样本（`POST /callset`，单样本文本也可直接提交）

紧凑表头：

```
#CHROM  POS  REF  ALT  S1    S2
chr1    100  A    T    0/1   1/1
```

VCF 风格表头（8 个固定列 + `FORMAT` + 样本列，基因型取自 `GT`）：

```
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1   S2
chr1   100 .  A   T   .    .      .    GT:DP 0/1:12 1/1:30
```

- 无表头时按首行样本列数推断，样本命名为 `sample_1 … sample_n`；之后列数不一致
  的行**不会被静默丢弃**，会报 `HEADER_SAMPLE_MISMATCH`。
- 表头样本名重复时报 `DUPLICATE_SAMPLE`，重复列自动改记为 `原名_2` 等，避免相互覆盖。

## 2. 归一化规则（等价写法收敛）

- 染色体别名归一：`1` == `chr1` == `CHR1` → `CHR1`。
- 等位基因大写、去空白；只接受非空 ACGT 碱基（N 不接受，直接拒绝）。
- 基因型归一到固定词表：`0/1`/`0|1`/`1/0`/`HET` → `0/1`；`1/1`/`1|1`/`HOM` → `1/1`；
  `0/0`/`REF`/`WT` → `0/0`；`.`/空/`UNCALLED` → `./.`。
- VCF 式左对齐：先去除共同末碱基，再逐列去除共同首碱基并把 POS 右移，每条等位
  基因始终保留至少一个锚定碱基。例如
  `chr1:100 ACG>AT` 与 `chr1:101 CG>T` 归一键相同：`CHR1:101:CG:T`。
- 变异身份键 `VariantKey = (chrom, pos, ref, alt)`；同键即等价写法。

## 3. 错误分类（可区分、带位置）

每条被拒绝的记录都返回 `reason` + `detail` + `location{block,line_no,raw}`，
`block` 形如 `sample:S1` 或 `multisample`，`line_no` 在该文本块内从 1 计数：

| reason | 触发条件 |
| --- | --- |
| `EMPTY_LINE` | 空内容（保留枚举，当前空行直接跳过） |
| `MISSING_FIELDS` | 列数不足，缺少 CHROM/POS/REF/ALT 等必需列 |
| `BAD_POSITION` | POS 不是正整数或 < 1 |
| `BAD_ALLELE` | REF/ALT 含非 ACGT 字符、为空，或 REF == ALT |
| `BAD_CHROM` | 染色体列为空 |
| `BAD_GENOTYPE` | 基因型不在可识别词表内（如 `9/9`） |
| `HEADER_SAMPLE_MISMATCH` | 数据行样本列数与表头（或首行推断）不一致 |
| `DUPLICATE_SAMPLE` | 表头中样本名重复 |
| `MALFORMED_LINE` | 列数过多、VCF `FORMAT` 缺 `GT` 等结构性损坏 |

HTTP 层面：畸形记录**不是**传输错误，返回 200 且 `valid=false`，逐条列出拒绝原因；
只有空提交/空样本名返回 422。

## 4. 合并与冲突规则

- 同一 `VariantKey` 的记录聚合为一个条目；携带样本去重，重复/等价写法**不会重复
  计数**（支持证据位置保留在 `support` 中）。
- 同一样本对同一变异的重复一致调用视为证据叠加，不算冲突。
- 输出按 `(染色体, POS, REF, ALT)` 确定性排序（1<2<…<10<X<Y<M<其它），重复提交结论
  与顺序稳定一致。
- 冲突**只标记、不用某一条覆盖另一条**，类型可区分：

| conflict type | 含义 |
| --- | --- |
| `GENOTYPE_CONFLICT` | 同一样本对同一变异同时出现 `0/1` 与 `1/1` |
| `REF_CONFLICT` | 同一样本在同一锚定位点同时被记为纯合 ALT 与纯合参考 `0/0` |
| `SAMPLE_CALL_CONFLICT` | 同一样本在同一位点对**不同 ALT** 都有 ALT 调用（各 ALT 都保留，各自标记） |

- 不同样本在同一位点携带不同 ALT 属正常情况，**不**判为冲突。
- `0/0` 不是 ALT 携带者（不出现在条目的 `samples` 中），但保留在
  `genotype_by_sample` 中用于冲突判定。

## 5. 并发、缓存与性能

- 每次提交的处理结果只存在于返回对象中，无样本级全局可变状态；并发处理有重叠的
  样本集合不会相互污染。
- 缓存键为 `sha256(sample + \\x00 + text)`：相同输入重复提交直接复用同一结论对象。
  缓存为有界线程安全 LRU（默认 128 条，`RLock` 保护，重活不持锁）。
- 解析/归约为 O(n log n)；仓库测试包含 1 万行输入的性能断言（本地远小于 5 秒）。

## 6. HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/`、`/health` | 服务信息 / 存活与 offline 标志 |
| POST | `/normalize` | body `{text, sample}`，单样本归一 |
| POST | `/callset` | body `{text}`，单/多样本归一+合并 |

响应字段：`content_hash`、`valid`、`accepted_count`、`entries[]`（含
`variant/chrom/pos/ref/alt/samples/genotype_by_sample/support/conflicts`）、
`rejected[]`、`conflicts[]`、`summary`。

启动：

```bash
uv run uvicorn main:app --host 127.0.0.1 --port 8000
# 或 uv run python main.py
```

## 7. 验证方法

```bash
uv sync
uv run pytest                      # 全部单测
uv run pytest -o log_cli=true      # 同时打印“输入摘要 / 判定依据”日志
```

测试覆盖（`tests/`）：

- 畸形输入：缺列、坏 POS/等位基因/基因型、样本列对不上、重名样本；断言错误分类
  可区分且带行号；
- 等价写法归一：染色体别名、VCF 左对齐、不同写法聚合成一个条目且不重复计数；
- 冲突记录：`GENOTYPE_CONFLICT`、`REF_CONFLICT`、`SAMPLE_CALL_CONFLICT`，以及
  “跨样本不同 ALT 不算冲突”；
- 并发与重复提交：16 线程 × 重叠样本集无污染、相同输入对象复用、乱序输入输出一致、
  LRU 有界淘汰、1 万行性能；
- 各单测通过 `logging` 打印输入摘要与判定依据（用 `-o log_cli=true` 查看）。
