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
- 缓存键为 `sha256(sample + \x00 + text)`：相同输入重复提交直接复用同一结论对象。
  缓存为有界线程安全 LRU（默认 128 条，`RLock` 保护，重活不持锁）。
- 每次提交（含增量，见第 6 节）还会生成不可变的**版本对象**存入独立的版本 LRU
  （默认 256 个），用 `version_id` 引用；版本读取只取短锁，不参与重计算。
- 解析/归约为 O(n log n)；仓库测试包含 1 万行输入的性能断言（本地远小于 5 秒）。

## 6. 增量更新与变更对比

样本是分批送达的，且可能补测/修正。无需重提全量样本：在上一轮返回的
`version_id` 基线上提交“只包含新增样本或待修正样本”的增量即可。

### 6.1 使用方式

`POST /increment`，请求体：

| 字段 | 说明 |
| --- | --- |
| `text` | 本次增量文本，**只带新增样本或需修正的样本**；同样支持单样本/多样本两种形态 |
| `baseline_version` | 上一轮（`/callset`、`/normalize` 或 `/increment`）返回的 `version_id` |
| `sample` | 可选。给出时按**单样本形态**解析（等价 `/normalize` 文本）；省略时按多样本形态解析（单样本 4 列文本也兼容） |

```bash
# 第一轮
curl -s localhost:8000/callset -H 'Content-Type: application/json' \
  -d '{"text":"#CHROM\tPOS\tREF\tALT\tS1\tS2\nchr1\t100\tA\tT\t0/1\t1/1\n"}'
# -> {"version_id":"v0...", ...}

# 增量：新增 S3 + 修正 S1（只给这两个样本即可）
curl -s localhost:8000/increment -H 'Content-Type: application/json' -d '{
  "baseline_version":"v0...",
  "text":"#CHROM\tPOS\tREF\tALT\tS3\tS1\nchr1\t100\tA\tT\t1/1\t1/1\nchr3\t300\tC\tG\t0/1\t0/1\n"}'

# 单样本形态的修正
curl -s localhost:8000/increment -H 'Content-Type: application/json' -d '{
  "baseline_version":"v1...", "sample":"S1", "text":"chr1 100 A T 0/1"}'
```

### 6.2 修正（整体替换）≠ 证据叠加

- **增量按样本整体替换**：本次携带的样本，先**删除**它此前的全部调用与可归属
  拒绝，再写入本次结论。因此同一样本重测后的新结论**不会**与旧结论被当成两份
  证据叠加（不会凭空产生 `GENOTYPE_CONFLICT`），该样本在各条目里的基因型、
  携带者身份、`support` 支持位置都以本次为准。
- 这与“同一次提交里重复一致调用算证据叠加、矛盾只标记不覆盖”的既有规则不冲突：
  若**本次**文本内部自相矛盾（如同时给 `0/1` 与 `1/1`），或修正样本的新结论与
  **其它样本**/它自己的其它新调用在同位点矛盾，仍照常按既有三类冲突标记，
  **不能用新结论直接盖过矛盾**。
- 本次增量没有碰到的旧条目、样本集合与排序**原样保留**（未触及条目直接复用，
  不重排、不消失）；样本的首见顺序稳定，新样本按表头顺序追加在后。
- 多样本形态里每个数据行都包含该批所有样本的单元格（含 `./.`），因此修正其中
  一个样本会一并带走它在该批各行的调用——这属于“以本次为准”的正确语义；只想
  改一个样本时用单样本形态（`sample` 字段）最精确。

### 6.3 版本标识与变更摘要

- `version_id`：对**最终整份 callset 内容**（条目、携带者与基因型、`support`
  溯源位置、拒绝、冲突、接受计数）做规范化 JSON 后取 `sha256`。它与请求形态
  解耦——**同一基线 + 同一批增量，无论提交几次，`version_id` 与 `change_summary`
  完全一致**；在同内容处派生出的版本还会复用同一个版本对象（幂等）。
- `parent_version_id`：上一轮基线 id，形成可追溯链条。
- `change_summary` 相对上一版本给出确定性、有序的差异：
  - `new_samples` / `corrected_samples` / `submitted_samples`：本次增量新增、
    修正、提交的样本；
  - `entries_added` / `entries_removed`：新增/消失的变异条目；
  - `entries_modified`：同一变异的变化，逐条给出 `samples_before/after` 与
    `sample_changes[]`，后者区分 `CARRIER_ADDED`（成为携带者）、
    `CARRIER_REMOVED`（不再携带）、`GENOTYPE_CHANGED`（如 `0/1`→`1/1`、与
    `0/0`/`./.` 之间的变化），并带 `before/after`；
  - `conflicts_appeared` / `conflicts_resolved`：新出现与已消除的冲突
    （类型/变异/样本/详情），冲突身份包含位置，确定性排序。

### 6.4 兼容范围与错误分类

- 两种输入形态、等价写法归一、全部既有 `reason`/冲突类型与排序规则对增量同样
  适用；列错位、坏等位基因、坏基因型、样本重名等仍是 **HTTP 200 + `valid=false`
  + 逐条 `rejected`（含 reason/line/sample）**，不会被当成正常增量悄悄吃入。
- 请求级、可区分、可复现的拒绝：
  - `404 BASELINE_NOT_FOUND`：`baseline_version` 不存在或已被 LRU 淘汰（版本库
    默认保留 256 个版本）；
  - `422 EMPTY_INCREMENT`：增量文本为空、样本名为空，或整批不含任何样本列
    （只有注释/表头）。
- 注意：增量里某样本的行**全部是坏数据**时，该样本本次没有有效结论，按“整体
  替换”它此前的有效调用会被移除、坏行进入 `rejected`——调用方据摘要与拒绝列表
  即可识别这是一次“修正失败”，而非被静默接受。

### 6.5 并发、缓存与性能

- 版本对象**不可变**且按内容寻址存于有界 LRU；解析/校验/重合并都在锁外进行
  （CPU 密集、仅用本地数据），所以大增量计算不会卡住对既有版本的读取
  （万行规模下基线读取为微秒级）。
- 同一基线并发提交多个增量（含样本重叠）时各自**派生独立子版本、互不污染、不丢
  样本**；因版本 id 是内容哈希，相同（基线, 增量）的并发请求确定地合流到同一
  版本对象，结果可复现。
- 增量采用**定向重合并**：只重算被修正/新增样本触及位点的条目与冲突，未触及条目
  直接复用；万行基线构建约 0.2s、单样本修正增量约 0.1s、万行新样本增量约 0.2s
  （本机数据，测试阈值更宽）。

## 7. HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/`、`/health` | 服务信息 / 存活与 offline 标志 |
| POST | `/normalize` | body `{text, sample}`，单样本归一 |
| POST | `/callset` | body `{text}`，单/多样本归一+合并 |
| POST | `/increment` | body `{text, baseline_version, sample?}`，在指定版本上增量更新 |

响应字段：`version_id`、`parent_version_id`、`change_summary`、`content_hash`、
`valid`、`accepted_count`、`entries[]`（含
`variant/chrom/pos/ref/alt/samples/genotype_by_sample/support/conflicts`）、
`rejected[]`（含可归属 `sample`）、`conflicts[]`、`summary`。

启动：

```bash
uv run uvicorn main:app --host 127.0.0.1 --port 8000
# 或 uv run python main.py
```

## 8. 验证方法

```bash
uv sync
uv run pytest                      # 全部单测
uv run pytest -o log_cli=true      # 同时打印“输入摘要 / 判定依据”日志
uv run pytest tests/test_increment.py tests/test_increment_service.py -o log_cli=true
```

测试覆盖（`tests/`）：

- 畸形输入：缺列、坏 POS/等位基因/基因型、样本列对不上、重名样本；断言错误分类
  可区分且带行号；
- 等价写法归一：染色体别名、VCF 左对齐、不同写法聚合成一个条目且不重复计数；
- 冲突记录：`GENOTYPE_CONFLICT`、`REF_CONFLICT`、`SAMPLE_CALL_CONFLICT`，以及
  “跨样本不同 ALT 不算冲突”；
- 并发与重复提交：16 线程 × 重叠样本集无污染、相同输入对象复用、乱序输入输出一致、
  LRU 有界淘汰、1 万行性能；
- **增量更新（本轮新增）**：
  - `test_increment.py`（纯逻辑层）：新增样本、修正整体替换而非证据叠加、未触及
    条目对象复用与顺序保持、条目/携带者增删、冲突出现与消除、同基线同增量可复现、
    畸形修正样本给出可区分原因且不静默生效、多样本列错位为批次级拒绝；
  - `test_increment_service.py`（服务层）：两种输入形态增量、重复增量返回同一版本
    对象的幂等性、`BASELINE_NOT_FOUND`/`EMPTY_INCREMENT` 可区分、40 个并发重叠
    增量不丢样本且版本序列可复现、万行规模耗时与基线读取不被阻塞；
  - `test_api.py` 增量用例：`/increment` 多样本/单样本形态 diff、HTTP 幂等、
    404/422 错误码与列错位 200/`valid=false` 响应；
- 各单测通过 `logging` 打印输入摘要与判定依据（用 `-o log_cli=true` 查看）。
