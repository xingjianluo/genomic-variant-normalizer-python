# 离线基因组变异记录处理服务（Genomic Variant Normalizer）

在原有 FastAPI 骨架上扩展而成。把多份样本的变异文本**归一、合并成一份可解释的
callset**。全过程只使用本地解析逻辑，**不访问任何外部服务 / 网络数据源**。

## 运行与验证

```bash
uv sync --group dev          # 安装依赖（含 pytest）
uv run pytest                 # 全部单测（日志会打印输入摘要与判定依据）
TEST_DEBUG=1 uv run pytest    # 更详细日志
uv run python main.py         # 启动服务（uvicorn，127.0.0.1 随机端口）
uv run uvicorn main:app --port 8000
```

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET  | `/health` | 健康检查，`offline: true` |
| POST | `/callsets` | 提交一批样本输入，返回**概览**（计数 + callset_id） |
| POST | `/callsets/inspect` | 提交并返回**完整明细**（变异、基因型汇总、冲突、拒绝记录） |
| GET  | `/callsets/{callset_id}` | 按 ID 取回已生成 callset（重复提交可直接复用） |

请求体为 JSON 数组，每个元素：

```json
{ "sample_id": "S1", "format": "single", "text": "chr1\t10\tA\tG\t0/1" }
```

* `format = "single"`：单样本形态；`sample_id` 必填，此时 `text` 每行
  `chrom <TAB> pos <TAB> ref <TAB> alt [<TAB> genotype]`（4 或 5 列）。
  若 `sample_id` 传空串，则每行首列作为样本名：
  `sample <TAB> chrom <TAB> pos <TAB> ref <TAB> alt [<TAB> genotype]`（5 或 6 列）。
* `format = "vcf"`：多样本 VCF 文本，以 `#CHROM ... FORMAT <samples...>` 表头声明
  样本列，每个样本的 GT 展开成一条 call；无样本列（8 列 sites-only）时整行归入
  样本 `UNKNOWN`、基因型 `UNKNOWN`。`##` 元信息行会被跳过。
* 空数组返回 `400`；畸形 JSON 返回 `422`。畸形**记录**不会让请求失败，而是进入
  响应的 `rejected_records`（HTTP 200）。

## 归一规则（等价写法收敛）

合并键为归一后的 `(chrom, pos, ref, alt)`：

* **染色体**：大小写/前缀归一，`1`、`CHR1` → `chr1`；`X/Y` 与 `M/MT` 等价
  （`chrM` → `chrMT`）；无法识别（如 `chrZ`）拒绝为 `INVALID_CHROM`。
* **等位基因**：一律大写，仅允许 `ACGTN`；为空、相同或符号等位（`<DEL>` 等）
  拒绝为 `INVALID_ALLELE`。
* **有界左移（无参考基因组 trimming）**：先去 REF/ALT 公共前缀（POS 随之左移，
  保留 1 个锚定碱基），再去公共后缀。可把带冗余锚定的等价写法收敛，例如
  `CA/CTA@100` 与 `A/TA@101` 归一键一致。跨位置的深层等价需要参考基因组，
  离线条件下不做推断（保守、可解释）。
* **基因型**：接受 `0/1`、`0|1`、`1/1`、`0/0`、`./.`、`.` 以及规范名
  `HOM_REF/HET/HOM_ALT/UNKNOWN`；多 ALT 位点按 GT 等位索引展开成多条 call，
  索引超出 ALT 数量拒绝为 `INVALID_GENOTYPE`；缺失/`./.` 记为 `UNKNOWN`。

## 错误分类（rejected_records，均可区分且带行号）

| reason | 触发场景 |
| --- | --- |
| `EMPTY_INPUT` | 输入为空或仅空白（`line=0`） |
| `MALFORMED_LINE` | VCF 数据行前缺少 `#CHROM` 表头；表头不是以 `#CHROM` 开头 |
| `WRONG_COLUMN_COUNT` | 行的列数不符合该形态要求；VCF 有 FORMAT 列却无样本值列 |
| `SAMPLE_COLUMN_MISMATCH` | VCF 样本值列数与表头声明的样本数不一致 |
| `DUPLICATE_SAMPLE_HEADER` | VCF 表头样本列重名 |
| `MISSING_FIELD` | chrom/pos/ref/alt 或样本名为空 |
| `INVALID_POSITION` | 位置不是整数或 `<1` |
| `INVALID_CHROM` | 染色体名无法识别 |
| `INVALID_ALLELE` | 等位基因含非法字符、为空、REF==ALT 或符号等位 |
| `INVALID_GENOTYPE` | 基因型无法解析，或 GT 等位索引超出 ALT 数量 |

每条拒绝包含 `source`（样本名或 `vcf`）、`line`（1 起始行号）、`reason`、
`detail`、`raw`（原始行），**绝不静默丢弃**。

## 合并与冲突规则

* **去重聚合**：同 `(样本, 变异键, 基因型)` 的重复记录（含等价写法收敛后的重复）
  只计一次；同一变异键下多个样本聚合成一个 callset 变异，不会重复计数。
* **GENOTYPE_CONFLICT（基因型冲突）**：同一样本在同一变异上出现互相矛盾的具体
  基因型——`HOM_REF` 与任意带 ALT 型别冲突；`HET` 与 `HOM_ALT` 冲突。
  `UNKNOWN` 不与任何型别冲突。**所有矛盾证据都保留在 `genotype_summary`，
  绝不任选一条覆盖另一条**。
* **SITE_CONFLICT（位点冲突）**：同一 `(chrom, pos)` 出现不同 REF 的记录。
  同 REF 不同 ALT 属于正常的多 ALT 位点，不算冲突；但位点 REF 存在歧义时，
  该位置的所有变异都会带 `SITE_CONFLICT` 标记与冲突 REF 列表。
* **稳定输出**：变异按 `(染色体规范顺序 1..22,X,Y,MT, pos, ref, alt)` 排序；
  拒绝记录按 `(source, line, reason)` 排序。同样输入无论提交顺序/次数，
  `callset_id`、输出顺序与结论完全一致。

## 并发、缓存与性能

* `callset_id` = 对排序后的输入集合 `(sample_id, format, text)` 取 SHA-256（16 位），
  因此与请求顺序无关，且重复提交天然命中同一结论。
* 每个 callset_id 一把异步锁：并发提交相同/重叠输入时只计算一次，其余请求复用
  缓存结果；计算过程是无共享可变状态的纯函数，样本集合重叠也不会互相污染。
* 解析/合并在线程池（默认 4 worker）中执行，不阻塞事件循环，限制并发重计算量；
  10k 行量级输入秒级完成（见 `tests/test_service.py` 的性能断言，上限 15s），
  同时并发的轻请求不会被饿死。

## 代码结构

```
variant_service/
  models.py        # 统一记录模型 / 枚举 / 错误原因
  normalization.py # 两种形态解析 + 染色体/等位/基因型归一 + 拒绝分类
  merge.py         # 去重聚合 + 基因型/位点冲突检测 + 稳定排序
  service.py       # 确定性 ID、缓存复用、按 key 锁、有界线程池
  api.py           # FastAPI 路由
tests/             # 畸形输入 / 等价归一 / 冲突 / 并发重复提交 测试
```
