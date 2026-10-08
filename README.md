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

## 5. 增量更新与版本对比

样本通常分批送达，同一样本还可能补测/修正。服务在不重放全量样本文本的前提下，
在**上一轮 callset** 上继续追加或替换样本，并产出新版本与相对上一版本的变更摘要。

### 5.1 版本标识（`version_id`）

- 每次提交（包括原 `/normalize`、`/callset`）都会得到一个 `version_id`：它是对
  **完整归一后记录状态**做 `sha256` 内容寻址得到的 64 位十六进制串（标签
  `callset-state-v1`）。
- 版本号**只对状态内容负责**，与提交顺序无关：两批互不重叠的样本，无论按什么顺序
  到达，最终 callset 的 `version_id` 完全一致。父版本 `parent_version` 只记录血缘，
  不参与版本号计算。
- 每个版本还带 `root_version`（血缘根）与 `sequence`（根为 0，每接受一次增量 +1）。
- 旧版本快照不可变，可用 `GET /version/{version_id}` 随时取回；版本存储为有界 LRU，
  **血缘头版本不会被淘汰**（它始终是合法基线）。

### 5.2 修正 ≠ 证据叠加（关键语义）

- 增量里出现的每个样本，其结论都会被**整批替换**：先移除该样本此前的全部归一记录与
  归属拒绝，再放入本次提交的全部记录。因此“样本修正”不会把新旧两条结论当成两份证据
  叠在一起（不会凭空制造 `GENOTYPE_CONFLICT`）。
- 同一次提交文本内部的重复/矛盾调用仍按既有规则处理：重复一致调用算证据、不重复计数；
  `0/1` 与 `1/1` 并存仍标 `GENOTYPE_CONFLICT`，同一锚定位点不同 ALT 仍标
  `SAMPLE_CALL_CONFLICT`。
- 修正样本如果与**其它样本**或该样本本次的其它调用矛盾，照样按第 4 节规则标记；
  新结论**不会**直接把旧冲突“盖没”——冲突是否消除，取决于重合并后的真实状态。
- 本次增量没有碰到的样本：其条目内容、携带者集合与相对排序保持原样（条目仍按
  `(染色体, POS, REF, ALT)` 全序输出）。

### 5.3 变更摘要 `changes`（相对所引用基线）

| 字段 | 含义 |
| --- | --- |
| `added_entries` / `removed_entries` | 新出现 / 消失的变异条目（`VariantKey` 列表） |
| `changed_entries[]` | 同一条目的变化：`samples_added`、`samples_removed`、`genotype_changes[]`（样本、旧基因型、新基因型） |
| `samples_added` / `samples_removed` | 在整个 callset 层面，新增/失去 ALT 携带身份的样本 |
| `new_conflicts` / `resolved_conflicts` | 相对基线新出现 / 已消除的冲突（完整冲突对象） |
| `new_rejections` / `resolved_rejections` | 新增 / 已消除的结构性拒绝（带原因与位置） |
| `identical` | 相对所引用基线是否无任何变化 |

同一基线 + 同一批增量，提交任意次，`version_id` 与 `changes` 都**逐字节一致**。

### 5.4 端点与输入形态兼容

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/increment/root/single` | body `{text, sample}`，建立版本 0（单样本形态） |
| POST | `/increment/root/multi` | body `{text}`，建立版本 0（紧凑/VCF/无表头形态均可） |
| POST | `/increment/single` | body `{text, sample, base_version}`，单样本增量/修正 |
| POST | `/increment/multi` | body `{text, base_version}`，多样本增量/修正 |
| GET | `/version/{version_id}` | 取回某个不可变版本 |

- 原 `/normalize`、`/callset` 的响应现在也带 `version_id`，可直接当作首轮基线引用；
  原有字段与错误分类完全不变。
- 单样本形态提交空（仅空白）文本，表示该样本**撤回全部调用**；HTTP 空文本仍返回 422。
- 畸形行（列错位、坏 POS/等位基因/基因型、重名样本等）沿用第 3 节可区分原因与行号，
  不会被当成正常增量吞掉；数据行级拒绝归属到该行样本，表头级拒绝归属到全部命名样本，
  样本修正后其历史拒绝随替换进入 `resolved_rejections`。

### 5.5 并发增量（乐观并发，确定结果，不丢样本）

针对同一基线并发提交多个增量时：

- 重建在锁外执行，**读取已有版本不会被大请求卡住**；提交时对血缘头做 CAS 校验。
- 样本**互不重叠**：自动重放到最新头版本并重试，内容寻址保证不同提交顺序收敛到同一
  最终版本，任何样本都不会丢。
- 样本**有重叠**且基线已被推进：返回 **HTTP 409**，`detail.reason =
  BASELINE_SUPERSEDED`，并给出 `current_version`（当前头版本）与
  `conflicting_samples`（冲突样本）——这是可区分、可复现的拒绝，客户端基于
  `current_version` 重新组织提交即可。
- 引用未知/已不存在的版本返回 **HTTP 404**，`reason = BASELINE_NOT_FOUND`；
  不含任何样本列的增量返回 422/`EMPTY_INCREMENT`。
- 重复提交（结果状态已存在）直接命中既有内容版本：对同一基线重复修正返回相同
  `version_id`，对当前头版本重复提交当前状态则返回 `identical=true` 的无操作结论。

### 5.6 调用示例

```bash
# 1) 首轮全量（老端点即可），记下 version_id
curl -s localhost:8000/callset -d '{"text":"#CHROM\tPOS\tREF\tALT\tS1\nchr1\t100\tA\tT\t0/1\n"}' \
  -H 'Content-Type: application/json'
# 2) 只送新增样本 S2，base_version 引用上一轮
curl -s localhost:8000/increment/single \
  -d '{"text":"chr2 200 G C 1/1","sample":"S2","base_version":"<上一轮 version_id>"}' \
  -H 'Content-Type: application/json'
# 3) 修正 S1（整批替换，不叠加）
curl -s localhost:8000/increment/single \
  -d '{"text":"chr1 100 A T 1/1","sample":"S1","base_version":"<第 2 步 version_id>"}' \
  -H 'Content-Type: application/json'
```

## 6. 并发、缓存与性能（原单次提交）

- 每次提交的处理结果只存在于返回对象中，无样本级全局可变状态；并发处理有重叠的
  样本集合不会相互污染。
- 缓存键为 `sha256(sample + \\x00 + text)`：相同输入重复提交直接复用同一结论。
  缓存为有界线程安全 LRU（默认 128 条，`RLock` 保护，重活不持锁）。
- 解析/归约为 O(n log n)；仓库测试包含 1 万行输入的性能断言（本地远小于 5 秒），
  增量场景另有“1 万行修正期间旧版本读取不被阻塞”的断言。

## 7. HTTP 接口汇总

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/`、`/health` | 服务信息 / 存活与 offline 标志 |
| POST | `/normalize` | body `{text, sample}`，单样本归一（响应带 `version_id`） |
| POST | `/callset` | body `{text}`，单/多样本归一+合并（响应带 `version_id`） |
| POST | `/increment/root/single`、`/increment/root/multi` | 建立版本化根 callset |
| POST | `/increment/single`、`/increment/multi` | 在 `base_version` 上应用增量 |
| GET | `/version/{version_id}` | 取回不可变版本快照 |

根/增量版本响应在第 3、4 节字段外，另含 `samples`、`touched_samples`、`sequence`、
`root_version`、`parent_version` 与 `changes`。

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
```

测试覆盖（`tests/`）：

- 畸形输入：缺列、坏 POS/等位基因/基因型、样本列对不上、重名样本；断言错误分类
  可区分且带行号；
- 等价写法归一：染色体别名、VCF 左对齐、不同写法聚合成一个条目且不重复计数；
- 冲突记录：`GENOTYPE_CONFLICT`、`REF_CONFLICT`、`SAMPLE_CALL_CONFLICT`，以及
  “跨样本不同 ALT 不算冲突”；
- 并发与重复提交：16 线程 × 重叠样本集无污染、相同输入结论复用、乱序输入输出一致、
  LRU 有界淘汰、1 万行性能；
- **增量更新**（`test_increment.py`、`test_increment_api.py`）：
  - 新增样本不扰动旧条目/排序、样本修正整批替换而非证据叠加、撤回样本使条目/携带者
    正确消失；
  - 变更摘要：新增/变化/移除条目、样本携带增减、基因型前后值、新冲突与已消除冲突、
    新增/已消除拒绝；
  - 版本可复现（同基线同增量多次提交版本号与摘要一致）、当前状态重复提交为
    `identical=true` 无操作；
  - 兼容两种输入形态；畸形增量仍返回可区分原因与位置，列错位不被吞掉；
  - 并发：64 个不相交样本并发增量不丢样本、且与顺序执行最终 `version_id` 完全一致；
    重叠样本并发增量要么命中确定内容版本，要么得到可区分、可复现的
    `BASELINE_SUPERSEDED`（409，含当前头版本与冲突样本）；未知基线为 404；
  - 1 万行修正重建期间，旧版本读取耗时远低于 0.1s、不被阻塞；
  - 各用例通过 `logging` 打印“判定依据[…]”日志（用 `-o log_cli=true` 查看）。
