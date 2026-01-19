# sim-reads 算法设计评估

> 版本: 0.8.16
> 日期: 2026-01-04

## 概述

`ecc sim-reads` 模拟从 eccDNA 样本通过 **Rolling Circle Amplification (RCA)** 扩增后进行测序的完整过程。本文档对其算法设计进行批判性评估。

---

## 核心流程

```
eccDNA 模板
    ↓
Step 1: Compartment 生成 (微滴/微孔分配)
    ↓
Step 2: RCA 分子图生成 (主干 + 分支结构)
    ↓
Step 3: Chimera 注入 (compartment 内分子互作)
    ↓
Step 4: 去支化/线性化 (文库制备)
    ↓
Step 4.5: Background DNA 生成 (可选)
    ↓
Step 5: 平台特异性 Read 生成 (NGS/HiFi/ONT)
    ↓
Step 6: 覆盖度补充 (确保每个 eccDNA 可检测)
    ↓
Step 7: 输出 (reads + truth + statistics)
```

---

## 优点

### 1. 生物学过程建模

**传统模拟器** 直接从序列随机切片生成 reads，忽略了扩增过程。

**本设计** 模拟了完整的分子生物学流程：

| 步骤 | 生物学意义 |
|------|-----------|
| Compartment | 微滴反应的空间隔离 |
| RCA 分子图 | φ29 聚合酶的滚环复制和链置换 |
| 分支结构 | 链置换产生的分叉 |
| Chimera | 同一反应空间内的模板跳转 |
| 去支化 | T7 核酸内切酶处理 |

生成的 reads 具有：
- 真实的重复结构（junction-spanning reads）
- 合理的 chimera 来源
- 分支来源的 reads

### 2. Compartment 模型

```python
P(chimera) = 1 - exp(-β × (N-1))
```

- 反映微滴/微孔反应的物理隔离
- Chimera 只在同一 compartment 内发生（符合真实情况）
- 避免了"全局随机 chimera"的不真实假设
- Compartment 大小直接影响 chimera 概率

### 3. 多平台差异化支持

| 平台 | 特点 | 模拟考虑 |
|------|------|----------|
| **NGS** | Paired-end 150bp | Insert size 分布, branch-chimera |
| **HiFi** | 长读长 ~15-20kb | 长度筛选 (≥5kb), CCS 错误模型 |
| **ONT** | 超长读长 | 全长分子测序, nanopore 错误谱 |

### 4. 覆盖度保证机制

采用分阶段生成策略：

```
Stage 1: 生成 50% 目标 reads
    ↓ 检查覆盖度
Stage 2: 补充覆盖不足的 eccDNA
    ↓
Stage 3: 生成剩余 reads
    ↓
PostFill: 最终检查和补充
```

确保：
- 每个 eccDNA 至少有 N× 覆盖（默认 5×）
- 低丰度 eccDNA 不会被遗漏
- 对 benchmark 测试至关重要

### 5. 完整的 Ground Truth

每条 read 包含完整的来源信息：

```
- source_ecc_ids: 来源 eccDNA 列表
- repeat_count: 重复次数
- has_chimera: 是否包含 inter-molecule chimera
- has_branch_chimera: 是否包含 branch chimera
- junction_covered_possible: 是否可能跨越 junction
- segments: 详细的序列来源坐标
```

便于下游工具的灵敏度/特异性评估。

### 6. 参数可配置

支持通过 YAML 配置文件调整：
- RCA 参数 (μ_R, σ_R, k_len)
- 分支参数 (B_rate, D_eff)
- Chimera 参数 (β, breakpoint 分布)
- 文库参数 (insert size, read length)
- 输出规模 (read count 或 coverage)

### 7. 并行化设计

- 多进程并行处理 RCA/Chimera/Debranch/Read 生成
- ecc_db 子集化减少数据传输（~95% 减少）
- 进度反馈机制

---

## 缺点 / 局限性

### 1. RCA 模型的简化

**当前模型**：
```python
ln(R) ~ N(μ_R - k_len × ln(L/L_ref), σ_R²)
```

**缺失的生物学因素**：

| 因素 | 真实影响 | 当前状态 |
|------|----------|----------|
| GC 含量 | 高 GC 区域聚合酶可能卡顿/脱落 | ❌ 未考虑 |
| 二级结构 | 发夹、G-quadruplex 阻碍复制 | ❌ 未考虑 |
| eccDNA 超螺旋 | 影响模板可及性和解旋 | ❌ 未考虑 |
| 引物结合效率 | 影响扩增起点和效率 | ❌ 假设均匀 |
| 序列 motif | 某些序列可能是聚合酶 pause site | ❌ 未考虑 |

### 2. 分支模型过于简单

**当前**：
```python
B ~ Poisson(B_rate × R)  # 分支数仅与重复次数相关
branch_length ~ LogNormal(μ_br, σ_br) × trunk_length
```

**问题**：
- 真实分支形成与链置换动力学、局部序列特征相关
- 分支长度分布缺乏实验数据支持
- 超分支（分支的分支）建模深度有限
- 分支锚点位置假设均匀分布，可能不真实

### 3. Chimera 断点分布

**当前**：
```python
breakpoint = uniform(0, trunk_length)  # 或 beta 分布
```

**问题**：
- 真实 chimera 断点可能偏好特定序列 motif
- 可能与 nick 位点、二级结构相关
- 缺乏真实数据的分布验证
- Beta 分布参数 (α, β) 是假设值

### 4. 文库制备模型粗糙

**去支化**：
```python
if random() < D_eff:  # 全局固定效率
    branch.debranched = True
```

**问题**：
- 真实 T7 核酸内切酶有序列偏好性
- 去支化效率可能与分支长度/结构相关
- Nick 到 break 的转化可能非随机
- Size selection 过程未精细建模

### 5. 错误模型局限

- 错误模型基于通用测序数据，非 eccDNA 特异性
- 缺乏 junction 附近的错误率变化建模
- 未考虑 PCR 扩增（如有）引入的错误
- HiFi CCS 模型可能需要更新

### 6. 性能瓶颈

```
大数据集处理仍然较慢：
- 复杂分子对象的创建和操作
- 序列字符串操作
- 跨进程数据序列化
```

即使经过并行优化，处理数万个 eccDNA 实例仍需较长时间。

### 7. 验证困难

**核心问题**：
- 真实 RCA 产物的分子结构难以直接观测
- 关键参数（μ_R, σ_R, B_rate, β 等）缺乏实验校准
- 只能通过下游结果间接验证
- 不同实验条件下参数可能差异很大

---

## 改进建议

### 短期改进（低成本高收益）

1. **GC 含量影响**
```python
def sample_repeat_count(self, ecc: EccDNA) -> int:
    gc = calculate_gc(ecc.seq)
    # GC 极端时降低扩增效率
    gc_factor = 1.0 - 0.5 * abs(gc - 0.45) ** 2
    base_R = self.sample_base_repeat(len(ecc.seq))
    return max(1, int(base_R * gc_factor))
```

2. **更详细的日志和统计**
   - 每个步骤的耗时
   - 分子结构统计（平均分支数、chimera 率等）

3. **序列化优化**
   - 只传递必要数据而非完整对象
   - 考虑使用共享内存

### 中期改进

1. **序列依赖的分支/chimera 模型**
   - 根据局部序列特征调整概率
   - 识别潜在的 pause site

2. **实验数据校准**
   - 收集真实 RCA 测序数据
   - 反推/拟合模型参数
   - 建立参数的置信区间

3. **文库制备精细化**
   - T7 切割偏好性
   - Size selection 模拟

### 长期改进

1. **结构感知模型**
   - 预测 eccDNA 二级结构（RNAfold 等）
   - 根据结构调整扩增效率和断点分布

2. **动力学模型**
   - 显式模拟 RCA 过程
   - 时间分辨的分支形成

3. **多条件支持**
   - 不同酶（φ29 vs Bst）
   - 不同反应条件（温度、时间、引物浓度）

---

## 适用场景

### 推荐使用

| 场景 | 说明 |
|------|------|
| 工具开发和调试 | 已知 ground truth，便于 debug |
| 灵敏度/特异性评估 | 可控制各种参数 |
| 参数探索 | 研究覆盖度、chimera 率等对检测的影响 |
| 方法比较 | 公平比较不同工具 |

### 谨慎使用

| 场景 | 注意事项 |
|------|----------|
| 定量分析验证 | 丰度关系可能与真实不完全一致 |
| 生物学结论 | 需结合真实数据验证 |

### 不推荐

| 场景 | 原因 |
|------|------|
| 替代真实实验 | 模型参数未经充分校准 |
| 作为唯一证据 | 模拟结果不能作为生物学发现的唯一依据 |

---

## 总体评价

| 维度 | 评分 | 说明 |
|------|------|------|
| 生物学合理性 | ⭐⭐⭐⭐☆ | 核心流程正确，细节有简化 |
| 工程实现 | ⭐⭐⭐☆☆ | 功能完整，性能有优化空间 |
| 可用性 | ⭐⭐⭐⭐☆ | 参数丰富，输出完整 |
| 可扩展性 | ⭐⭐⭐⭐☆ | 模块化设计，易于添加新功能 |
| 可验证性 | ⭐⭐☆☆☆ | 缺乏实验校准数据 |

**总结**：

该设计在概念层面是合理的，捕捉了 RCA 扩增的核心特征。对于 **benchmark 测试和工具开发** 是有价值的工具。但由于缺乏实验数据校准，不能完全替代真实数据，建议与真实实验数据结合使用。

---

## 参考

- φ29 DNA Polymerase: Blanco & Salas, 1996
- Rolling Circle Amplification: Fire & Xu, 1995
- eccDNA detection methods: Various, 2020-2024
