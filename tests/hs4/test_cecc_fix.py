#!/usr/bin/env python3
"""
测试 cecc_build 的修复方案：正确处理重叠比对

核心问题：cecc_build 在遇到重叠比对时直接丢弃整个 query，
而不是选择最佳的非重叠比对子集。

修复方案：使用贪心算法选择非重叠的最佳比对记录
"""

import pandas as pd
import numpy as np
from pathlib import Path


def select_non_overlapping_alignments(group: pd.DataFrame,
                                       overlap_tolerance: int = 10) -> pd.DataFrame:
    """
    从重叠的比对记录中选择最佳的非重叠子集。

    使用贪心算法：按 alignment_length 降序排列，
    逐个选择不与已选记录重叠的比对。

    Args:
        group: 单个 query 的所有比对记录
        overlap_tolerance: 允许的重叠容忍度 (bp)

    Returns:
        选择的非重叠比对记录
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
            # 计算重叠
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


def preprocess_query_alignments(df: pd.DataFrame,
                                 overlap_tolerance: int = 10) -> pd.DataFrame:
    """
    预处理所有 query 的比对记录，选择非重叠子集。

    Args:
        df: 包含所有比对记录的 DataFrame
        overlap_tolerance: 允许的重叠容忍度 (bp)

    Returns:
        处理后的 DataFrame
    """
    processed_parts = []

    for query_id, group in df.groupby('query_id', sort=False):
        processed = select_non_overlapping_alignments(group, overlap_tolerance)
        processed_parts.append(processed)

    if not processed_parts:
        return df.iloc[0:0].copy()

    result = pd.concat(processed_parts, ignore_index=True)
    return result


def find_circular_fixed(group_df: pd.DataFrame,
                        edge_tol: int = 20,
                        min_coverage: float = 0.9) -> dict:
    """
    修复后的环形路径检测。

    Args:
        group_df: 单个 query 的比对记录（已预处理）
        edge_tol: 允许的边缘间隙 (bp)
        min_coverage: 最小覆盖率要求

    Returns:
        检测结果 dict，或 None
    """
    if len(group_df) < 2:
        return None

    # 按 q_start 排序
    sorted_df = group_df.sort_values('q_start').reset_index(drop=True)

    # 获取序列长度（双倍化后的长度约为 2 * length）
    cons_len = float(sorted_df['length'].iloc[0])
    doubled_len = cons_len * 2

    # 检查连续性并构建路径
    path = [0]
    total_coverage = float(sorted_df.iloc[0]['alignment_length'])

    for idx in range(1, len(sorted_df)):
        prev = sorted_df.iloc[path[-1]]
        cur = sorted_df.iloc[idx]

        gap = float(cur['q_start']) - float(prev['q_end'])

        # 允许一定的间隙
        if gap > edge_tol:
            # 间隙过大，停止
            break

        # 负间隙说明有重叠（预处理后不应该发生大的重叠）
        if gap < -edge_tol:
            continue  # 跳过这个重叠记录

        path.append(idx)
        total_coverage += float(cur['alignment_length'])

    # 计算覆盖率
    coverage_ratio = total_coverage / doubled_len if doubled_len > 0 else 0

    if coverage_ratio < min_coverage:
        return None

    return {
        'path': path,
        'coverage': coverage_ratio,
        'total_aligned': total_coverage,
        'num_segments': len(path),
    }


def test_single_query(query_records: pd.DataFrame, query_id: str):
    """测试单个 query 的修复效果"""
    print(f"\n{'='*60}")
    print(f"测试 {query_id}")
    print(f"{'='*60}")

    print(f"\n原始比对记录数: {len(query_records)}")

    # 预处理：选择非重叠子集
    processed = select_non_overlapping_alignments(query_records, overlap_tolerance=10)
    print(f"预处理后记录数: {len(processed)}")

    # 显示处理前后的对比
    print(f"\n预处理后的比对记录 (按 q_start 排序):")
    sorted_proc = processed.sort_values('q_start')
    for i, (_, row) in enumerate(sorted_proc.iterrows()):
        print(f"  {i}: q={int(row['q_start'])}-{int(row['q_end'])} "
              f"({int(row['alignment_length'])}bp) -> {row['chr']}")

    # 计算间隙
    print(f"\n间隙分析:")
    prev_end = 0
    max_gap = 0
    for _, row in sorted_proc.iterrows():
        gap = row['q_start'] - prev_end
        if prev_end > 0:
            print(f"  gap = {int(gap)} bp")
            max_gap = max(max_gap, abs(gap))
        prev_end = row['q_end']

    print(f"\n最大间隙: {max_gap} bp")

    # 尝试检测环形
    result = find_circular_fixed(processed, edge_tol=50)

    if result:
        print(f"\n✓ 检测成功!")
        print(f"  路径长度: {result['num_segments']} 个片段")
        print(f"  覆盖率: {result['coverage']*100:.1f}%")
        print(f"  总比对长度: {result['total_aligned']:.0f} bp")
    else:
        print(f"\n✗ 检测失败")

    return result


def main():
    """主测试函数"""
    # 加载测试数据
    data_dir = Path('/Users/yaoxinzhang/Documents/eccToolkit/tests/hs4')
    unclass_file = data_dir / 'um_classify.unclassified.csv'
    cecc_build_file = data_dir / 'cecc_build.csv'
    truth_file = data_dir / 'hs1_sim.chimeric.bed'

    print("加载数据...")
    unclass_df = pd.read_csv(unclass_file)
    cecc_build_df = pd.read_csv(cecc_build_file)

    # 加载真值
    truth_df = pd.read_csv(truth_file, sep='\t', comment='#',
                           names=['chrom','start','end','name','length','strand','type','fragments'],
                           header=None)
    if truth_df.iloc[0]['chrom'] == 'chrom':
        truth_df = truth_df.iloc[1:].reset_index(drop=True)

    # 提取 CeccDNA 相关记录
    def extract_truth_id(query_id):
        parts = str(query_id).split('.')
        return parts[0] if len(parts) >= 1 else query_id

    unclass_df['truth_id'] = unclass_df['query_id'].apply(extract_truth_id)
    cecc_build_df['truth_id'] = cecc_build_df['query_id'].apply(extract_truth_id)

    cecc_unclass = unclass_df[unclass_df['truth_id'].str.startswith('CeccDNA')]

    # 找出漏检的 CeccDNA
    detected_truth_ids = set(cecc_build_df[cecc_build_df['truth_id'].str.startswith('CeccDNA')]['truth_id'])
    all_cecc_truth_ids = set(truth_df['name'])
    missed_truth_ids = all_cecc_truth_ids - detected_truth_ids

    print(f"\n原始检测结果:")
    print(f"  CeccDNA 真值总数: {len(all_cecc_truth_ids)}")
    print(f"  被检测到: {len(detected_truth_ids)}")
    print(f"  漏检: {len(missed_truth_ids)}")

    # 测试几个漏检的样本
    print("\n" + "="*60)
    print("测试漏检样本的修复效果")
    print("="*60)

    # 选几个漏检的 query 测试
    test_samples = ['CeccDNA_000001', 'CeccDNA_000008', 'CeccDNA_000009']

    fixed_count = 0
    for truth_id in test_samples:
        if truth_id not in missed_truth_ids:
            continue

        # 获取这个 truth_id 的所有 query
        sample_records = cecc_unclass[cecc_unclass['truth_id'] == truth_id]

        # 选第一个 query 测试
        first_qid = sample_records['query_id'].iloc[0]
        query_records = sample_records[sample_records['query_id'] == first_qid]

        result = test_single_query(query_records, first_qid)
        if result:
            fixed_count += 1

    # 批量测试所有漏检的 CeccDNA
    print("\n" + "="*60)
    print("批量测试所有漏检的 CeccDNA")
    print("="*60)

    fixed_total = 0
    still_missed = 0

    for truth_id in missed_truth_ids:
        sample_records = cecc_unclass[cecc_unclass['truth_id'] == truth_id]
        if sample_records.empty:
            continue

        # 测试所有 query，只要有一个成功就算修复
        query_fixed = False
        for qid in sample_records['query_id'].unique():
            query_records = sample_records[sample_records['query_id'] == qid]
            processed = select_non_overlapping_alignments(query_records, overlap_tolerance=10)

            if len(processed) >= 2:
                result = find_circular_fixed(processed, edge_tol=50)
                if result and result['coverage'] >= 0.9:
                    query_fixed = True
                    break

        if query_fixed:
            fixed_total += 1
        else:
            still_missed += 1

    print(f"\n修复结果统计:")
    print(f"  原漏检数: {len(missed_truth_ids)}")
    print(f"  修复后可检测: {fixed_total}")
    print(f"  仍然漏检: {still_missed}")
    print(f"  修复率: {fixed_total/len(missed_truth_ids)*100:.1f}%")

    # 计算新的召回率
    new_detected = len(detected_truth_ids) + fixed_total
    new_recall = new_detected / len(all_cecc_truth_ids) * 100
    print(f"\n新召回率: {new_detected}/{len(all_cecc_truth_ids)} = {new_recall:.1f}%")
    print(f"(原召回率: {len(detected_truth_ids)/len(all_cecc_truth_ids)*100:.1f}%)")


if __name__ == '__main__':
    main()
