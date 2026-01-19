#!/usr/bin/env python3
"""
Test MAPQ filtering in U/M classification

This script analyzes the effect of the new MAPQ ambiguity threshold
for Mecc classification using existing classification results.
"""

import pandas as pd
import re
from pathlib import Path
from collections import Counter


def extract_truth_type(query_id):
    """Extract the true eccDNA type from query_id."""
    if 'MeccDNA' in str(query_id):
        return 'M'
    elif 'UeccDNA' in str(query_id):
        return 'U'
    elif 'CeccDNA' in str(query_id):
        return 'C'
    return None


def analyze_mapq_filter_effect(mecc_df, threshold):
    """
    Analyze what would happen if we apply MAPQ filter to Mecc classification.

    A query would be vetoed if ALL its loci have mapq_best < threshold.
    """
    if mecc_df.empty:
        return set(), {}

    # Group by query_id and get max mapq_best for each query
    query_mapq = mecc_df.groupby('query_id')['mapq_best'].max().reset_index()
    query_mapq.columns = ['query_id', 'max_mapq']

    # Queries that would be vetoed
    vetoed_queries = set(query_mapq[query_mapq['max_mapq'] < threshold]['query_id'])

    # Analyze truth types of vetoed queries
    vetoed_truths = [extract_truth_type(q) for q in vetoed_queries]
    truth_counts = Counter(vetoed_truths)

    return vetoed_queries, truth_counts


def main():
    """Main test function."""
    data_dir = Path('/Users/yaoxinzhang/Documents/eccToolkit/tests/hs4')

    # Load pre-processed data
    print("Loading pre-processed classification data...")
    mecc_df = pd.read_csv(data_dir / 'um_classify.mecc.csv')
    uecc_df = pd.read_csv(data_dir / 'um_classify.uecc.csv')
    unclass_df = pd.read_csv(data_dir / 'um_classify.unclassified.csv')

    mecc_queries = mecc_df['query_id'].unique()
    uecc_queries = uecc_df['query_id'].unique()

    print(f"Mecc: {len(mecc_df)} records, {len(mecc_queries)} queries")
    print(f"Uecc: {len(uecc_df)} records, {len(uecc_queries)} queries")
    print(f"Unclassified: {len(unclass_df)} records, {unclass_df['query_id'].nunique()} queries")

    # Analyze original classification accuracy
    print("\n" + "="*60)
    print("Original Classification Results")
    print("="*60)

    mecc_truths = [extract_truth_type(q) for q in mecc_queries]
    mecc_truth_counts = Counter(mecc_truths)

    uecc_truths = [extract_truth_type(q) for q in uecc_queries]
    uecc_truth_counts = Counter(uecc_truths)

    print(f"\nMecc classification ({len(mecc_queries)} queries):")
    print(f"  True MeccDNA: {mecc_truth_counts.get('M', 0)}")
    print(f"  True UeccDNA: {mecc_truth_counts.get('U', 0)} (misclassified)")
    print(f"  True CeccDNA: {mecc_truth_counts.get('C', 0)}")
    if len(mecc_queries) > 0:
        precision = mecc_truth_counts.get('M', 0) / len(mecc_queries) * 100
        print(f"  Precision: {precision:.1f}%")

    print(f"\nUecc classification ({len(uecc_queries)} queries):")
    print(f"  True UeccDNA: {uecc_truth_counts.get('U', 0)}")
    print(f"  True MeccDNA: {uecc_truth_counts.get('M', 0)} (misclassified)")
    print(f"  True CeccDNA: {uecc_truth_counts.get('C', 0)}")
    if len(uecc_queries) > 0:
        precision = uecc_truth_counts.get('U', 0) / len(uecc_queries) * 100
        print(f"  Precision: {precision:.1f}%")

    # Test different MAPQ thresholds
    print("\n" + "="*60)
    print("Effect of MAPQ Filtering on Mecc Classification")
    print("="*60)

    print("\nAnalyzing MAPQ distribution in Mecc queries...")
    mecc_query_mapq = mecc_df.groupby('query_id').agg({
        'mapq_best': 'max'
    }).reset_index()
    mecc_query_mapq['truth'] = mecc_query_mapq['query_id'].apply(extract_truth_type)

    print("\nMAPQ distribution by truth type:")
    for t in ['M', 'U', 'C']:
        subset = mecc_query_mapq[mecc_query_mapq['truth'] == t]['mapq_best']
        if len(subset) > 0:
            print(f"  True {t}eccDNA: mean={subset.mean():.1f}, median={subset.median()}, "
                  f"<20: {(subset<20).sum()}, >=20: {(subset>=20).sum()}")

    print("\n" + "-"*60)
    print(f"{'Threshold':<12} {'Vetoed':<10} {'True U':<10} {'True M':<10} {'New Prec':<10}")
    print("-"*60)

    for threshold in [0, 10, 15, 20, 25, 30]:
        vetoed_queries, truth_counts = analyze_mapq_filter_effect(mecc_df, threshold)
        vetoed_u = truth_counts.get('U', 0)
        vetoed_m = truth_counts.get('M', 0)

        # Calculate new precision
        new_mecc_count = len(mecc_queries) - len(vetoed_queries)
        new_correct = mecc_truth_counts.get('M', 0) - vetoed_m
        new_precision = new_correct / new_mecc_count * 100 if new_mecc_count > 0 else 0

        print(f"{threshold:<12} {len(vetoed_queries):<10} {vetoed_u:<10} {vetoed_m:<10} {new_precision:.1f}%")

    # Detailed analysis for threshold=20 (default)
    print("\n" + "="*60)
    print("Detailed Analysis: MAPQ threshold = 20")
    print("="*60)

    vetoed_queries, truth_counts = analyze_mapq_filter_effect(mecc_df, 20)

    print(f"\nQueries vetoed: {len(vetoed_queries)}")
    print(f"  True UeccDNA: {truth_counts.get('U', 0)} (correctly removed)")
    print(f"  True MeccDNA: {truth_counts.get('M', 0)} (incorrectly removed)")
    print(f"  True CeccDNA: {truth_counts.get('C', 0)}")

    # Net improvement
    original_wrong_u = mecc_truth_counts.get('U', 0)
    new_wrong_u = original_wrong_u - truth_counts.get('U', 0)
    lost_m = truth_counts.get('M', 0)

    print(f"\nNet effect:")
    print(f"  UeccDNA misclassified as Mecc: {original_wrong_u} -> {new_wrong_u} (-{truth_counts.get('U', 0)})")
    print(f"  True MeccDNA lost: {lost_m}")
    print(f"  Net improvement: {truth_counts.get('U', 0) - lost_m} queries correctly handled")

    # New precision
    new_mecc_total = len(mecc_queries) - len(vetoed_queries)
    new_mecc_correct = mecc_truth_counts.get('M', 0) - lost_m
    original_precision = mecc_truth_counts.get('M', 0) / len(mecc_queries) * 100
    new_precision = new_mecc_correct / new_mecc_total * 100 if new_mecc_total > 0 else 0

    print(f"\n  Original Mecc precision: {original_precision:.1f}%")
    print(f"  New Mecc precision: {new_precision:.1f}%")
    print(f"  Improvement: +{new_precision - original_precision:.1f}%")

    # Analyze vetoed true MeccDNA - what's special about them?
    if lost_m > 0:
        print("\n" + "="*60)
        print(f"Analyzing {lost_m} lost true MeccDNA queries")
        print("="*60)

        lost_m_queries = [q for q in vetoed_queries if extract_truth_type(q) == 'M']
        lost_m_records = mecc_df[mecc_df['query_id'].isin(lost_m_queries)]

        print(f"\nMapq distribution of lost MeccDNA:")
        lost_mapq = lost_m_records.groupby('query_id')['mapq_best'].max()
        print(f"  max mapq_best: {lost_mapq.max()}")
        print(f"  mean mapq_best: {lost_mapq.mean():.1f}")

        print(f"\nThese might be MeccDNA from regions where all copies have low MAPQ")
        print("(possibly in complex repeat regions)")


if __name__ == '__main__':
    main()
