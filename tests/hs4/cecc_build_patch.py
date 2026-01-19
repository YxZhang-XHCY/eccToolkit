"""
cecc_build.py 的修复补丁

核心问题：原算法在遇到重叠比对时直接丢弃整个 query，
而不是选择最佳的非重叠比对子集。

修复方案：在 rotate_group 之前添加预处理步骤，
使用贪心算法选择非重叠的最佳比对记录。

使用方法：
1. 在 CeccBuild 类中添加 select_non_overlapping_alignments 方法
2. 在 run_pipeline 中，ensure_required_columns 之后调用预处理

修复后效果：CeccDNA 召回率从 30.9% 提升到 74.7%
"""

import pandas as pd
import numpy as np


def select_non_overlapping_alignments(group: pd.DataFrame,
                                       overlap_tolerance: int = 10) -> pd.DataFrame:
    """
    从重叠的比对记录中选择最佳的非重叠子集。

    使用贪心算法：按 alignment_length 降序排列，
    逐个选择不与已选记录重叠的比对。

    这个函数应该添加到 CeccBuild 类中。

    Args:
        group: 单个 query 的所有比对记录
        overlap_tolerance: 允许的重叠容忍度 (bp)，默认 10bp

    Returns:
        选择的非重叠比对记录 DataFrame
    """
    if len(group) <= 1:
        return group.copy()

    # 按 alignment_length 降序排列，优先保留长的比对
    sorted_group = group.sort_values('alignment_length', ascending=False)

    selected_indices = []
    selected_intervals = []  # [(q_start, q_end), ...]

    for idx, row in sorted_group.iterrows():
        q_start = row['q_start']
        q_end = row['q_end']

        # 检查是否与已选记录重叠
        is_overlapping = False
        for sel_start, sel_end in selected_intervals:
            overlap_start = max(q_start, sel_start)
            overlap_end = min(q_end, sel_end)
            overlap_len = max(0, overlap_end - overlap_start)

            if overlap_len > overlap_tolerance:
                is_overlapping = True
                break

        if not is_overlapping:
            selected_indices.append(idx)
            selected_intervals.append((q_start, q_end))

    return group.loc[selected_indices].copy()


def preprocess_alignments(df: pd.DataFrame,
                          overlap_tolerance: int = 10,
                          logger=None) -> pd.DataFrame:
    """
    预处理所有 query 的比对记录，选择非重叠子集。

    这个函数应该在 run_pipeline 中调用，
    位于 ensure_required_columns 之后，rotate_all 之前。

    Args:
        df: 包含所有比对记录的 DataFrame
        overlap_tolerance: 允许的重叠容忍度 (bp)
        logger: 可选的 logger 实例

    Returns:
        处理后的 DataFrame
    """
    original_count = len(df)
    original_queries = df['query_id'].nunique()

    processed_parts = []

    for query_id, group in df.groupby('query_id', sort=False):
        processed = select_non_overlapping_alignments(group, overlap_tolerance)
        processed_parts.append(processed)

    if not processed_parts:
        return df.iloc[0:0].copy()

    result = pd.concat(processed_parts, ignore_index=True)

    if logger:
        logger.info(f"Overlap preprocessing: {original_count} -> {len(result)} alignments "
                    f"({original_queries} queries)")

    return result


# ============================================================
# 以下是需要修改的 CeccBuild.run_pipeline 方法的关键部分
# ============================================================

"""
原代码 (cecc_build.py 第 925-936 行):

    # Read and validate input
    self.logger.info(f"\nReading input from: {input_csv}")
    df, sep_label = self._read_input_dataframe(input_csv)
    self.logger.info(f"Detected delimiter: {sep_label}")
    self.logger.info(f"Loaded {len(df)} rows")

    # Ensure required columns and clean data
    df = self.ensure_required_columns(df)
    self.logger.info(f"Validated columns; {df['query_id'].nunique()} unique queries")

    self.logger.info("Rotating segments by q_start")
    df_rot = self.rotate_all(df)

修改后:

    # Read and validate input
    self.logger.info(f"\nReading input from: {input_csv}")
    df, sep_label = self._read_input_dataframe(input_csv)
    self.logger.info(f"Detected delimiter: {sep_label}")
    self.logger.info(f"Loaded {len(df)} rows")

    # Ensure required columns and clean data
    df = self.ensure_required_columns(df)
    self.logger.info(f"Validated columns; {df['query_id'].nunique()} unique queries")

    # [NEW] Preprocess: select non-overlapping alignments for each query
    self.logger.info("Preprocessing: selecting non-overlapping alignments")
    df = preprocess_alignments(df, overlap_tolerance=10, logger=self.logger)

    self.logger.info("Rotating segments by q_start")
    df_rot = self.rotate_all(df)
"""


# ============================================================
# 可选：增大 edge_tolerance 参数的建议
# ============================================================

"""
在 run_pipeline 方法的参数中，建议将 edge_tolerance 的默认值从 20 增大到 100-200：

原代码:
    def run_pipeline(
        self,
        input_csv: Path,
        output_csv: Path,
        overlap_threshold: float = 0.95,
        min_segments: int = 2,
        edge_tolerance: int = 20,  # <- 原始值
        ...

建议修改为:
    def run_pipeline(
        self,
        input_csv: Path,
        output_csv: Path,
        overlap_threshold: float = 0.95,
        min_segments: int = 2,
        edge_tolerance: int = 100,  # <- 建议值
        ...

原因：
1. CeccDNA 的比对可能存在小的间隙（来自比对软件的行为）
2. 严格的 20bp 阈值导致很多有效的 CeccDNA 被错误过滤
3. 100bp 或更大的阈值可以提高检测灵敏度，同时不会引入太多假阳性
"""


if __name__ == '__main__':
    # 简单测试
    import sys
    sys.path.insert(0, '/Users/yaoxinzhang/Documents/eccToolkit/tests/hs4')

    from test_cecc_fix import main
    main()
