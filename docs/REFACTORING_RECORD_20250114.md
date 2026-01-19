# eccToolkit 重构记录 (2025-01-14)

## 背景问题

### 1. Lite 模式 ONT 无法生成
**原因**：代码使用 PBSIM2 参数格式（`--hmm_model`），但未检测工具是否支持。

**解决**：添加 PBSIM2 兼容性检测，不支持时跳过 ONT 并给出警告。

### 2. Full 模式与 Lite 模式架构冗余
**原问题**：
- Lite 模式：简化 RCA + 外部工具 (ART/PBSIM) 生成 reads
- Full 模式：复杂 RCA + **内置 Python library generator** 生成 reads

**核心洞察**：两者的本质差异只在 RCA 复杂度，reads 生成完全可以统一。

---

## 重构方案

### 新架构

```
Full 模式:
  eccDNA FASTA
      ↓
  Step 1: Compartment 分室
      ↓
  Step 2: RCA 分子生成 (分支/超分支/动力学)
      ↓
  Step 3: Chimera 注入
      ↓
  Step 4: Debranch 去分支
      ↓
  Step 4.5: 背景 DNA 生成 (可选)
      ↓
  导出 CSV (含 eccDNA 来源追踪)
      ↓
  fqsim (外部工具: ART/PBSIM)
      ↓
  NGS/HiFi/ONT reads

Lite 模式:
  eccDNA FASTA
      ↓
  简化 RCA (固定长度扩增)
      ↓
  导出 CSV
      ↓
  fqsim (外部工具)
      ↓
  reads
```

**核心改动**：Full 模式彻底移除内置 library generator，统一使用 fqsim。

---

## 具体修改

### 1. `src/ecctoolkit/simulate/lite.py`

**添加 PBSIM2 检测** (第 17-20 行):
```python
def _check_tool_available(tool_name: str) -> bool:
    """检测外部工具是否可用"""
    return shutil.which(tool_name) is not None
```

**修改 `sim_fastq_ont()` 方法** (第 760-790 行):
- 在生成 ONT reads 前检测 PBSIM2 是否支持 `--hmm_model` 参数
- 不支持则警告并跳过

### 2. `src/ecctoolkit/simulate/full/rca_readsim/debranch.py`

**添加 `to_dataframe()` 和 `to_csv()` 方法** (第 698-769 行):

```python
def to_dataframe(self, include_coverage: bool = True, include_truth: bool = True) -> pd.DataFrame:
    """将分子池转换为 DataFrame 格式，供 fqsim 使用"""
    # 输出字段：
    # - id: 分子 ID
    # - seq: 序列
    # - length: 长度
    # - tempcov: 覆盖度
    # - ecc_ids: eccDNA 来源 ID (用 | 分隔)
    # - repeat_count: RCA 重复次数
    # - source_ecc_length: 来源 eccDNA 长度
    # - has_chimera: 是否有嵌合
    # - is_background: 是否是背景 DNA
    # - source_graph_id: RCA 图 ID
    # - background_region: 背景 DNA 区域

def to_csv(self, output_path: str, ...) -> str:
    """将分子池导出为 CSV 文件"""
```

### 3. `src/ecctoolkit/simulate/full/reads.py`

**彻底重写**：
- 原文件：1268 行
- 新文件：537 行
- **删除**：内置 library generator 代码 (~700 行)
- **保留**：RCA pipeline (Step 1-4.5)
- **新增**：调用 fqsim 生成 reads

关键代码 (第 485-524 行):
```python
# ========== Generate reads using external tools (ART/PBSIM) ==========
logger.info("Step 5: Generating reads using external tools (ART/PBSIM)...")

from ..lite import fqsim

# Export molecule pool to CSV
pool_csv_path = output_path / f"{prefix}pool.lib.csv"
merged_pool.to_csv(str(pool_csv_path), include_coverage=True, include_truth=True)

# Call fqsim to generate reads
fqsim(
    sample=sample or "sample",
    csv=str(pool_csv_path),
    path=str(output_path),
    ...
)
```

### 4. `src/ecctoolkit/simulate/pipeline.py`

**修改 `_run_full_readsim()`**:
- 添加 `lite = self.config.readsim.lite` (从 lite 配置获取外部工具参数)
- 传递外部工具参数给 `run_read_simulation()`

### 5. `src/ecctoolkit/simulate/unified_config.py`

- 移除 `use_external_tools` 字段 (不再需要，默认就是外部工具)

### 6. `src/ecctoolkit/cli.py`

- 移除 `--use-external-tools` 选项

---

## 测试结果

### 测试命令
```bash
python -m src.ecctoolkit.cli simulate \
    -r demo_region_test/demo_region_test.all.fa \
    -o _test_full_external \
    --skip-region \
    --input-eccdna demo_region_test/demo_region_test.all.fa \
    -u 10 -t 2 --seed 42 \
    --readsim-mode full \
    --skip-ont --skip-hifi
```

### 测试输出
```
Step 1: Generating compartments...
Step 2: Generating RCA molecules...
Step 3: Injecting chimeras...
Step 4: Debranching...
Step 4.5: Generating background linear DNA...
  Molecule pool: 11380 eccDNA + 1328 background
Step 5: Generating reads using external tools (ART/PBSIM)...
  Exported molecule pool to: sim_ecc_pool.lib.csv
  Total molecules: 12708
  Calling fqsim (external tools)...
  External tools read generation completed.
```

### 生成文件
| 文件 | 大小 | 说明 |
|------|------|------|
| `sim_ecc_pool.lib.csv` | 4.4GB | 分子池（含 eccDNA 追踪） |
| `sim_ecc.NGS.R1.fastq` | 422MB | NGS R1 reads |
| `sim_ecc.NGS.R2.fastq` | 422MB | NGS R2 reads |

### CSV 格式示例
```
id                      seq         length  tempcov  ecc_ids       repeat_count  has_chimera  is_background
mol_0_1_trunk_frag0     GGCGGACT... 5230    12.5     ecc_001       15            False        False
mol_0_1_trunk_frag1     TTGCCATG... 8100    8.3      ecc_001|002   20            True         False
bg_mol_0                AATCGTAG... 3200    5.0                    1             False        True
```

---

## 优势

1. **速度提升**：外部工具 (ART/PBSIM) 比纯 Python 实现快得多
2. **效果可靠**：ART/PBSIM 是成熟的模拟器，误差模型更准确
3. **代码简化**：删除 700+ 行冗余代码，维护成本降低
4. **eccDNA 追踪**：CSV 完整记录每个分子的来源，支持下游分析
5. **架构统一**：Full/Lite 的差异只在 RCA，reads 生成统一

---

## 待办事项

1. [ ] 可以考虑删除不再使用的模块：
   - `src/ecctoolkit/simulate/full/rca_readsim/library/` (内置 library generator)
   - `src/ecctoolkit/simulate/full/rca_readsim/error_models/` (内置误差模型)

2. [ ] 进一步简化：考虑是否合并 Lite 和 Full 为单一 pipeline，通过参数控制 RCA 复杂度

3. [ ] 测试 HiFi 和 ONT 生成是否正常

---

## 相关文件

- 本文档：`docs/REFACTORING_RECORD_20250114.md`
- 之前的优化计划：`docs/FULL_MODE_OPTIMIZATION.md`
- 之前的计划文件：`~/.claude/plans/purring-sniffing-puzzle.md`
