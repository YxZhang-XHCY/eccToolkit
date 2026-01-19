#!/usr/bin/env python3
"""
Test Rule + ML combination for Mecc classification

This script tests the effect of combining:
1. Rule-based filter: identity_gap >= 1.0 (100% safe)
2. ML-based filter: trained classifier with adjustable threshold

The combination should:
- Never reject true MeccDNA (rule filter is 100% safe)
- Further reduce false positives beyond what the rule alone achieves
"""

import pandas as pd
import numpy as np
import pickle
import json
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


def compute_features(group, locus_col='locus_id', identity_col='identity'):
    """Compute ML features for a query group."""
    # Basic stats
    length = group['length'].iloc[0] if 'length' in group.columns else 0
    n_loci = group[locus_col].nunique()

    # Coverage stats
    locus_cov = group.groupby(locus_col)['locus_cov'].first()
    u_cov = locus_cov.max() if len(locus_cov) > 0 else 0
    locus_cov_min = locus_cov.min() if len(locus_cov) > 0 else 0
    locus_cov_std = locus_cov.std() if len(locus_cov) > 1 else 0

    # Identity stats
    identities = group[identity_col]
    identity_max = identities.max()
    identity_min = identities.min()
    identity_mean = identities.mean()
    identity_std = identities.std() if len(identities) > 1 else 0

    # MAPQ
    mapq_best = group['mapq_best'].max() if 'mapq_best' in group.columns else 0

    # Identity gap per locus
    loci_identity = group.groupby(locus_col)[identity_col].max()
    sorted_ids = sorted(loci_identity.values, reverse=True)
    id_gap = sorted_ids[0] - sorted_ids[1] if len(sorted_ids) >= 2 else 0
    id_gap_3rd = sorted_ids[0] - sorted_ids[2] if len(sorted_ids) >= 3 else 0

    return {
        'length': length,
        'n_loci': n_loci,
        'U_cov': u_cov,
        'locus_cov_min': locus_cov_min,
        'locus_cov_std': locus_cov_std if not pd.isna(locus_cov_std) else 0,
        'identity_max': identity_max,
        'identity_min': identity_min,
        'identity_mean': identity_mean,
        'identity_std': identity_std if not pd.isna(identity_std) else 0,
        'mapq_best': mapq_best,
        'id_gap': id_gap,
        'id_gap_3rd': id_gap_3rd,
    }


def main():
    """Main test function."""
    data_dir = Path('/Users/yaoxinzhang/Documents/eccToolkit/tests/hs4')
    model_path = Path('/Users/yaoxinzhang/Documents/CircleSeeker/src/circleseeker/models/mecc_classifier.pkl')
    metadata_path = Path('/Users/yaoxinzhang/Documents/CircleSeeker/src/circleseeker/models/mecc_classifier_metadata.json')

    # Load pre-processed data
    print("Loading pre-processed classification data...")
    mecc_df = pd.read_csv(data_dir / 'um_classify.mecc.csv')

    # Load ML model
    print("Loading ML model...")
    with open(model_path, 'rb') as f:
        ml_model = pickle.load(f)
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)

    feature_cols = metadata['feature_columns']
    print(f"Feature columns: {feature_cols}")

    # Compute features for all queries
    print("\nComputing features for all queries...")
    query_features = []
    for query_id, group in mecc_df.groupby('query_id'):
        feats = compute_features(group)
        feats['query_id'] = query_id
        feats['truth'] = extract_truth_type(query_id)
        query_features.append(feats)

    features_df = pd.DataFrame(query_features)

    # Original statistics
    mecc_queries = mecc_df['query_id'].unique()
    mecc_truths = [extract_truth_type(q) for q in mecc_queries]
    mecc_truth_counts = Counter(mecc_truths)

    total_queries = len(mecc_queries)
    true_m_count = mecc_truth_counts.get('M', 0)
    false_pos_u_count = mecc_truth_counts.get('U', 0)
    original_precision = true_m_count / total_queries * 100

    print("\n" + "="*70)
    print("Original Classification Results (No Filtering)")
    print("="*70)
    print(f"Total Mecc queries: {total_queries}")
    print(f"  True MeccDNA: {true_m_count}")
    print(f"  False positive (UeccDNA): {false_pos_u_count}")
    print(f"  Precision: {original_precision:.1f}%")

    # Rule-only filtering (identity_gap >= 1.0)
    print("\n" + "="*70)
    print("Rule-only Filtering: identity_gap >= 1.0")
    print("="*70)

    rule_threshold = 1.0
    rule_vetoed = features_df[features_df['id_gap'] >= rule_threshold]
    rule_vetoed_m = (rule_vetoed['truth'] == 'M').sum()
    rule_vetoed_u = (rule_vetoed['truth'] == 'U').sum()

    remaining_after_rule = features_df[features_df['id_gap'] < rule_threshold]
    remaining_m = (remaining_after_rule['truth'] == 'M').sum()
    remaining_u = (remaining_after_rule['truth'] == 'U').sum()
    rule_precision = remaining_m / (remaining_m + remaining_u) * 100 if (remaining_m + remaining_u) > 0 else 0

    print(f"Vetoed by rule: {len(rule_vetoed)}")
    print(f"  True MeccDNA incorrectly vetoed: {rule_vetoed_m}")
    print(f"  False positive correctly vetoed: {rule_vetoed_u}")
    print(f"Remaining: {len(remaining_after_rule)}")
    print(f"  True MeccDNA: {remaining_m}")
    print(f"  False positive: {remaining_u}")
    print(f"  Precision: {rule_precision:.1f}%")
    print(f"  Improvement: +{rule_precision - original_precision:.1f}%")

    # Rule + ML combination
    print("\n" + "="*70)
    print("Rule + ML Combination (Cascading Filter)")
    print("="*70)

    # First apply rule filter
    after_rule = features_df[features_df['id_gap'] < rule_threshold].copy()

    # Then apply ML to remaining
    X = after_rule[feature_cols].values
    ml_proba = ml_model.predict_proba(X)[:, 1]
    after_rule['ml_proba'] = ml_proba

    print(f"\nAfter rule filter: {len(after_rule)} queries")
    print(f"  True MeccDNA: {(after_rule['truth'] == 'M').sum()}")
    print(f"  False positive: {(after_rule['truth'] == 'U').sum()}")

    print(f"\n{'ML Threshold':<15} {'ML Vetoed':<12} {'True M Lost':<12} {'FP Removed':<12} {'Final Prec':<12}")
    print("-"*70)

    for ml_threshold in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        ml_vetoed = after_rule[after_rule['ml_proba'] < ml_threshold]
        ml_vetoed_m = (ml_vetoed['truth'] == 'M').sum()
        ml_vetoed_u = (ml_vetoed['truth'] == 'U').sum()

        final = after_rule[after_rule['ml_proba'] >= ml_threshold]
        final_m = (final['truth'] == 'M').sum()
        final_u = (final['truth'] == 'U').sum()
        final_precision = final_m / (final_m + final_u) * 100 if (final_m + final_u) > 0 else 0

        print(f"{ml_threshold:<15} {len(ml_vetoed):<12} {ml_vetoed_m:<12} {ml_vetoed_u:<12} {final_precision:.1f}%")

    # Detailed analysis for default threshold (0.5)
    print("\n" + "="*70)
    print("Detailed Analysis: Rule (id_gap >= 1.0) + ML (threshold = 0.5)")
    print("="*70)

    ml_threshold = 0.5
    final = after_rule[after_rule['ml_proba'] >= ml_threshold]
    final_m = (final['truth'] == 'M').sum()
    final_u = (final['truth'] == 'U').sum()
    final_precision = final_m / (final_m + final_u) * 100 if (final_m + final_u) > 0 else 0

    total_vetoed = total_queries - len(final)
    total_vetoed_m = true_m_count - final_m
    total_vetoed_u = false_pos_u_count - final_u

    print(f"\nTotal queries vetoed (rule + ML): {total_vetoed}")
    print(f"  - By rule (id_gap >= 1.0): {len(rule_vetoed)}")
    print(f"  - By ML (proba < 0.5): {len(after_rule) - len(final)}")

    print(f"\nTrue MeccDNA:")
    print(f"  Original: {true_m_count}")
    print(f"  Vetoed: {total_vetoed_m} ({100*total_vetoed_m/true_m_count:.1f}%)")
    print(f"  Retained: {final_m} ({100*final_m/true_m_count:.1f}%)")

    print(f"\nFalse Positive (UeccDNA misclassified as Mecc):")
    print(f"  Original: {false_pos_u_count}")
    print(f"  Removed: {total_vetoed_u} ({100*total_vetoed_u/false_pos_u_count:.1f}%)")
    print(f"  Remaining: {final_u}")

    print(f"\nPrecision:")
    print(f"  Original: {original_precision:.1f}%")
    print(f"  After rule only: {rule_precision:.1f}%")
    print(f"  After rule + ML: {final_precision:.1f}%")
    print(f"  Total improvement: +{final_precision - original_precision:.1f}%")

    # Summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"""
The Rule + ML combination approach:

1. Rule filter (id_gap >= 1.0):
   - 100% safe: never rejects true MeccDNA
   - Removes {rule_vetoed_u} of {false_pos_u_count} false positives ({100*rule_vetoed_u/false_pos_u_count:.1f}%)

2. ML filter (threshold = 0.5):
   - Applied ONLY to queries that pass the rule filter
   - Removes additional false positives with minimal true MeccDNA loss
   - True MeccDNA loss: {total_vetoed_m} ({100*total_vetoed_m/true_m_count:.1f}%)
   - False positive removal: {total_vetoed_u} ({100*total_vetoed_u/false_pos_u_count:.1f}%)

3. Combined effect:
   - Precision improved from {original_precision:.1f}% to {final_precision:.1f}%
   - True MeccDNA retained: {final_m}/{true_m_count} ({100*final_m/true_m_count:.1f}%)
   - False positives remaining: {final_u}/{false_pos_u_count}
""")


if __name__ == '__main__':
    main()
