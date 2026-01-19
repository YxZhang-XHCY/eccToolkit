# eccToolkit Full 模式性能优化方案

## 概述

本文档描述 Full 模式读段模拟的性能瓶颈分析和优化方案。

**目标**: 50-70% 总体加速
**技术栈**: Python + NumPy
**典型数据规模**: 10K-100K eccDNA

---

## 一、当前架构

```
Full 模式 Pipeline:
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
  导出 CSV (含 eccDNA 来源追踪)  ← 瓶颈 #1
      ↓
  fqsim (外部工具: ART/PBSIM)    ← 瓶颈 #2
      ↓
  NGS/HiFi/ONT reads
```

**关键改动**（2025-01-14 重构）：Full 模式移除内置 library generator，统一使用 fqsim。

---

## 二、性能瓶颈分析

### 2.1 Profiling 结果（2025-01-14）

测试条件：10 eccDNA, 100x coverage, NGS only

| 阶段 | 时间 | 占比 | 备注 |
|------|------|------|------|
| 配置加载 | ~9 秒 | ~5% | 包含配置解析、参数初始化 |
| Step 1-4 (RCA) | ~0.3 秒 | **<1%** | **非常快！** |
| Step 4.5 (背景 DNA) | ~2 秒 | ~1% | |
| **CSV 导出** | **~55 秒** | **~30%** | **瓶颈 #1** |
| **fqsim (ART)** | **预估 20+ 分钟** | **~65%** | **瓶颈 #2** |

### 2.2 主要瓶颈

#### 瓶颈 #1: CSV 导出（`to_csv()`）

**文件**: `src/ecctoolkit/simulate/full/rca_readsim/debranch.py:698-769`

**问题**:
1. 每个分子调用 `get_sequence()` 重建完整序列
2. 内部循环遍历 segments，调用 `get_circular_substr()` 和 `reverse_complement()`
3. 字符串拼接低效
4. 生成 4.1GB CSV 文件（12,708 分子，平均 330KB/分子）

**调用链**:
```
to_csv()
  → to_dataframe()
    → get_sequence(mol)  # 每个分子
      → mol.get_sequence(ecc_db)
        → ecc.get_circular_substr()  # 每个 segment
        → reverse_complement()       # 如果是反向链
```

#### 瓶颈 #2: fqsim 逐个调用 ART

**文件**: `src/ecctoolkit/simulate/lite.py` (fqsim 函数)

**问题**:
1. 每个分子单独调用一次 `art_illumina`
2. 12,708 个分子 = 12,708 次进程启动
3. 进程启动开销累积巨大

---

## 三、优化方案

### Phase 1: CSV 导出优化（预期 20-30% 总体提升）

#### 1.1 NumPy 向量化反向互补

**文件**: `src/ecctoolkit/simulate/full/rca_readsim/seq_utils.py`

**当前实现**:
```python
def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', ...}
    return "".join(complement.get(base, 'N') for base in reversed(seq))
```

**优化实现**:
```python
import numpy as np

# 预计算查找表（模块级常量）
_COMPLEMENT_TABLE = np.zeros(256, dtype=np.uint8)
_COMPLEMENT_TABLE[ord('A')] = ord('T')
_COMPLEMENT_TABLE[ord('T')] = ord('A')
_COMPLEMENT_TABLE[ord('G')] = ord('C')
_COMPLEMENT_TABLE[ord('C')] = ord('G')
_COMPLEMENT_TABLE[ord('a')] = ord('t')
_COMPLEMENT_TABLE[ord('t')] = ord('a')
_COMPLEMENT_TABLE[ord('g')] = ord('c')
_COMPLEMENT_TABLE[ord('c')] = ord('g')
_COMPLEMENT_TABLE[ord('N')] = ord('N')
_COMPLEMENT_TABLE[ord('n')] = ord('n')

def reverse_complement_fast(seq: str) -> str:
    """NumPy 向量化反向互补"""
    seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    rc_bytes = _COMPLEMENT_TABLE[seq_bytes[::-1]]
    return rc_bytes.tobytes().decode('ascii')
```

**预期提升**: 3-5%

---

#### 1.2 EccDNA NumPy 序列缓存

**文件**: `src/ecctoolkit/simulate/full/rca_readsim/models.py`

```python
@dataclass
class EccDNA:
    id: str
    seq: str
    weight: float = 1.0

    # 懒加载 NumPy 表示（不参与序列化）
    _seq_bytes: Optional[np.ndarray] = field(default=None, repr=False, compare=False)
    _seq_bytes_rc: Optional[np.ndarray] = field(default=None, repr=False, compare=False)

    @property
    def seq_bytes(self) -> np.ndarray:
        """懒加载 NumPy 字节数组"""
        if self._seq_bytes is None:
            self._seq_bytes = np.frombuffer(self.seq.encode('ascii'), dtype=np.uint8)
        return self._seq_bytes

    @property
    def seq_bytes_rc(self) -> np.ndarray:
        """预计算反向互补字节数组"""
        if self._seq_bytes_rc is None:
            from .seq_utils import _COMPLEMENT_TABLE
            self._seq_bytes_rc = _COMPLEMENT_TABLE[self.seq_bytes[::-1]]
        return self._seq_bytes_rc

    def get_circular_substr_fast(self, offset: int, length: int) -> str:
        """NumPy 加速的环状序列提取"""
        L = len(self.seq)
        offset = offset % L

        if length <= L:
            if offset + length <= L:
                return self.seq_bytes[offset:offset + length].tobytes().decode('ascii')
            else:
                part1 = self.seq_bytes[offset:]
                part2 = self.seq_bytes[:length - len(part1)]
                return np.concatenate([part1, part2]).tobytes().decode('ascii')
        else:
            cycles_needed = (offset + length + L - 1) // L
            extended = np.tile(self.seq_bytes, cycles_needed)
            return extended[offset:offset + length].tobytes().decode('ascii')
```

**预期提升**: 10-15%

---

#### 1.3 批量序列构建

**文件**: `src/ecctoolkit/simulate/full/rca_readsim/debranch.py`

**问题**: 当前 `to_dataframe()` 逐个分子构建序列

**优化**: 预分配 + 批量处理

```python
def to_dataframe_fast(self, include_coverage: bool = True, include_truth: bool = True) -> 'pd.DataFrame':
    """批量构建 DataFrame（优化版）"""
    import pandas as pd

    n = len(self.molecules)

    # 预分配列表
    ids = [None] * n
    seqs = [None] * n
    lengths = [None] * n

    # 批量处理
    for idx, mol in enumerate(self.molecules):
        ids[idx] = mol.molecule_id
        seqs[idx] = self.get_sequence(mol)  # 可进一步优化
        lengths[idx] = len(seqs[idx])

    data = {'id': ids, 'seq': seqs, 'length': lengths}

    if include_coverage:
        data['tempcov'] = [float(self.weights[i] * n) for i in range(n)]

    # ... truth 字段类似处理

    return pd.DataFrame(data)
```

**预期提升**: 5-10%

---

### Phase 2: fqsim 批处理优化（预期 30-40% 总体提升）

#### 2.1 合并分子 FASTA

**当前问题**: 每个分子单独一个 FASTA 文件，单独调用 ART

**优化方案**: 将多个分子合并为批次

```python
def fqsim_batch(molecules, batch_size=1000, ...):
    """批量处理分子"""
    for batch_start in range(0, len(molecules), batch_size):
        batch = molecules[batch_start:batch_start + batch_size]

        # 合并为单个 FASTA
        batch_fasta = create_batch_fasta(batch)

        # 单次 ART 调用
        run_art(batch_fasta, ...)

        # 拆分输出
        split_reads(batch, ...)
```

**预期提升**: 30-40%（减少进程启动开销）

---

#### 2.2 并行 ART 调用

**当前问题**: 串行处理

**优化方案**: 多进程并行

```python
from concurrent.futures import ProcessPoolExecutor

def fqsim_parallel(molecules, num_workers=4, ...):
    """并行调用 ART"""
    batches = split_into_batches(molecules, num_workers)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_art_batch, batch) for batch in batches]
        results = [f.result() for f in futures]

    merge_results(results)
```

**预期提升**: 额外 10-20%（取决于 CPU 核心数）

---

### Phase 3: 架构优化（长期）

#### 3.1 考虑不在 CSV 中存储完整序列

**问题**: 4.1GB CSV 主要是序列数据

**方案 A**: 只存储元数据，fqsim 按需提取序列
```
CSV: id, length, tempcov, segments_json
fqsim: 根据 segments_json 从 ecc_db 提取序列
```

**方案 B**: 使用二进制格式（如 pickle/parquet）存储分子池

---

#### 3.2 共享内存优化

**文件**: `src/ecctoolkit/simulate/full/rca_readsim/shared_pool.py` (新建)

用于并行场景下避免序列化 ecc_db。

---

## 四、实施计划

### Phase 1 任务清单（优先）

| 任务 | 文件 | 优先级 | 状态 |
|------|------|--------|------|
| 添加 `_COMPLEMENT_TABLE` 常量 | `seq_utils.py` | P0 | [ ] |
| 添加 `reverse_complement_fast()` | `seq_utils.py` | P0 | [ ] |
| EccDNA 添加 `seq_bytes` 属性 | `models.py` | P0 | [ ] |
| EccDNA 添加 `get_circular_substr_fast()` | `models.py` | P0 | [ ] |
| 更新 `get_sequence()` 使用快速方法 | `models.py` | P1 | [ ] |
| 优化 `to_dataframe()` | `debranch.py` | P1 | [ ] |

### Phase 2 任务清单

| 任务 | 文件 | 优先级 | 状态 |
|------|------|--------|------|
| fqsim 批量合并 FASTA | `lite.py` | P0 | [ ] |
| fqsim 并行 ART 调用 | `lite.py` | P1 | [ ] |
| 减少临时文件 I/O | `lite.py` | P2 | [ ] |

---

## 五、预期收益

| 阶段 | 优化项 | 预期提升 |
|------|--------|----------|
| Phase 1 | NumPy 序列操作 + 批量构建 | **20-30%** |
| Phase 2 | fqsim 批处理 + 并行 | **30-40%** |
| **总计** | | **50-70%** |

---

## 六、测试与验证

### 6.1 基准测试命令

```bash
# 小规模测试
python -m cProfile -s cumulative -o profile.pstats -m src.ecctoolkit.cli simulate \
    -r demo_region_test/demo_region_test.all.fa \
    -o /tmp/benchmark \
    --skip-region \
    --input-eccdna demo_region_test/demo_region_test.all.fa \
    -u 100 -t 1 --seed 42 \
    --readsim-mode full \
    --skip-ont --skip-hifi

# 查看结果
python -c "import pstats; p = pstats.Stats('profile.pstats'); p.sort_stats('cumulative').print_stats(30)"
```

### 6.2 正确性验证

优化后需验证：
1. 输出 FASTQ 序列正确
2. eccDNA 来源追踪正确
3. 与优化前结果一致（相同 seed）

---

## 更新日志

- 2025-01-14 v2: 基于新架构重新分析瓶颈，更新优化方案
  - 移除已删除的 library/error_models 相关内容
  - 新增 CSV 导出和 fqsim 优化方案
- 2025-01-14 v1: 初始版本（已过时）
