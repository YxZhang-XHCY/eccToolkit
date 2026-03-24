# Full RCA 模式 (开发中)

这个目录包含完整的 RCA (滚环扩增) 模拟模型代码。

## 状态

**当前状态**: 开发中 (WIP - Work In Progress)

Full 模式试图精确模拟 RCA 的生物学过程，包括：
- 支链 (branch) 结构生成
- 嵌合体 (chimera) 注入
- 去支化 (debranch) 处理
- 分子打断 (fragmentation)

## 已知问题

在代码审阅中发现了以下问题需要修复：

### 严重问题
1. NGS/HiFi 分子池缺少 chimera 信息传递
2. HiFi 打断会丢失 <5kb 的小 eccDNA
3. 覆盖度抽样逻辑可能导致偏差

### 中等问题
4. NGS 打断点使用均匀分布，不够真实
5. Chimera 断点重叠处理不够健壮
6. 背景 DNA 数据库传递可能有问题

## 目录结构

```
_wip_full/
├── __init__.py
├── bed_to_fasta.py      # BED 转 FASTA
├── reads.py             # 主入口
└── rca_readsim/         # RCA 模拟核心
    ├── models.py        # 数据结构
    ├── config.py        # 配置
    ├── rca_engine.py    # RCA 引擎
    ├── chimera.py       # 嵌合体注入
    ├── debranch.py      # 去支化
    ├── fragmentation.py # 打断
    ├── kinetics.py      # 动力学模型
    └── ...
```

## 恢复使用

当 Full 模式修复完成后，可以通过以下步骤恢复：

1. 将 `_wip_full` 重命名为 `full`
2. 更新 `simulate/__init__.py` 添加导入
3. 更新 `simulate/pipeline.py` 恢复 Full 模式处理
4. 更新 `simulate/unified_config.py` 恢复默认配置

## 相关文档

- Lite 模式在 `simulate/lite.py` 中实现
- 统一配置在 `simulate/unified_config.py` 中定义
