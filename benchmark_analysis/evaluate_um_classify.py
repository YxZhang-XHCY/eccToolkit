#!/usr/bin/env python3
"""
Evaluate UeccDNA/MeccDNA classification performance against ground truth.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


def load_ground_truth(truth_file: Path) -> dict[str, str]:
    """Load ground truth and extract true type for each read."""
    df = pd.read_csv(truth_file, sep='\t')

    read_to_type = {}
    for _, row in df.iterrows():
        read_id = row['read_id']
        # Extract type from source_mol_id (UeccDNA_xxx, MeccDNA_xxx, CeccDNA_xxx)
        source_mol = row['source_mol_id']
        if source_mol.startswith('UeccDNA'):
            true_type = 'U'
        elif source_mol.startswith('MeccDNA'):
            true_type = 'M'
        elif source_mol.startswith('CeccDNA'):
            true_type = 'C'
        else:
            true_type = 'X'  # Unknown
        read_to_type[read_id] = true_type

    return read_to_type


def load_classification_results(uecc_file: Path, mecc_file: Path, unclass_file: Path) -> dict[str, str]:
    """Load classification results and extract predicted type for each read."""
    read_to_pred = {}

    # Load Uecc
    if uecc_file.exists() and uecc_file.stat().st_size > 0:
        df_uecc = pd.read_csv(uecc_file)
        if not df_uecc.empty and 'reads' in df_uecc.columns:
            for reads in df_uecc['reads'].dropna():
                # reads column contains the original read ID
                read_to_pred[str(reads)] = 'U'

    # Load Mecc
    if mecc_file.exists() and mecc_file.stat().st_size > 0:
        df_mecc = pd.read_csv(mecc_file)
        if not df_mecc.empty and 'reads' in df_mecc.columns:
            for reads in df_mecc['reads'].dropna():
                read_to_pred[str(reads)] = 'M'

    # Load Unclassified
    if unclass_file.exists() and unclass_file.stat().st_size > 0:
        df_unclass = pd.read_csv(unclass_file)
        if not df_unclass.empty and 'reads' in df_unclass.columns:
            for reads in df_unclass['reads'].dropna().unique():
                if str(reads) not in read_to_pred:
                    read_to_pred[str(reads)] = 'X'  # Unclassified

    return read_to_pred


def compute_metrics(true_labels: list, pred_labels: list, positive_class: str):
    """Compute precision, recall, F1 for a specific class."""
    tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == positive_class and p == positive_class)
    fp = sum(1 for t, p in zip(true_labels, pred_labels) if t != positive_class and p == positive_class)
    fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == positive_class and p != positive_class)
    tn = sum(1 for t, p in zip(true_labels, pred_labels) if t != positive_class and p != positive_class)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
        'Precision': precision, 'Recall': recall, 'F1': f1
    }


def build_confusion_matrix(true_labels: list, pred_labels: list, classes: list) -> pd.DataFrame:
    """Build confusion matrix."""
    matrix = defaultdict(lambda: defaultdict(int))

    for t, p in zip(true_labels, pred_labels):
        matrix[t][p] += 1

    # Convert to DataFrame
    data = []
    for true_class in classes:
        row = [matrix[true_class][pred_class] for pred_class in classes]
        data.append(row)

    df = pd.DataFrame(data, index=[f'True_{c}' for c in classes], columns=[f'Pred_{c}' for c in classes])
    return df


def analyze_errors(read_to_type: dict, read_to_pred: dict):
    """Analyze specific error cases."""
    errors = {
        'U_as_M': [],  # UeccDNA misclassified as MeccDNA
        'M_as_U': [],  # MeccDNA misclassified as UeccDNA
        'U_as_X': [],  # UeccDNA unclassified
        'M_as_X': [],  # MeccDNA unclassified
        'C_as_U': [],  # CeccDNA misclassified as UeccDNA
        'C_as_M': [],  # CeccDNA misclassified as MeccDNA
    }

    for read_id, true_type in read_to_type.items():
        pred_type = read_to_pred.get(read_id, 'X')

        if true_type == 'U' and pred_type == 'M':
            errors['U_as_M'].append(read_id)
        elif true_type == 'M' and pred_type == 'U':
            errors['M_as_U'].append(read_id)
        elif true_type == 'U' and pred_type == 'X':
            errors['U_as_X'].append(read_id)
        elif true_type == 'M' and pred_type == 'X':
            errors['M_as_X'].append(read_id)
        elif true_type == 'C' and pred_type == 'U':
            errors['C_as_U'].append(read_id)
        elif true_type == 'C' and pred_type == 'M':
            errors['C_as_M'].append(read_id)

    return errors


def main():
    base_dir = Path('/Users/yaoxinzhang/Documents/eccToolkit/benchmark_analysis')

    # Load ground truth
    truth_file = base_dir / 'Demo_1400.HiFi.truth.tsv'
    read_to_type = load_ground_truth(truth_file)

    # Load classification results
    uecc_file = base_dir / 'um_classify.uecc.csv'
    mecc_file = base_dir / 'um_classify.mecc.csv'
    unclass_file = base_dir / 'um_classify.unclassified.csv'
    read_to_pred = load_classification_results(uecc_file, mecc_file, unclass_file)

    # Summary statistics
    print("=" * 70)
    print("Ground Truth Distribution:")
    type_counts = defaultdict(int)
    for t in read_to_type.values():
        type_counts[t] += 1
    for t in ['U', 'M', 'C']:
        print(f"  {t}eccDNA: {type_counts[t]:,} reads")
    print(f"  Total: {len(read_to_type):,} reads")

    print("\n" + "=" * 70)
    print("Classification Results Distribution:")
    pred_counts = defaultdict(int)
    for p in read_to_pred.values():
        pred_counts[p] += 1
    for p in ['U', 'M', 'X']:
        print(f"  Predicted {p}: {pred_counts[p]:,} reads")
    print(f"  Total classified: {len(read_to_pred):,} reads")

    # Build aligned lists (only consider reads present in both)
    common_reads = set(read_to_type.keys()) & set(read_to_pred.keys())
    print(f"\n  Reads in common: {len(common_reads):,}")

    true_labels = [read_to_type[r] for r in common_reads]
    pred_labels = [read_to_pred[r] for r in common_reads]

    # Confusion matrix
    print("\n" + "=" * 70)
    print("Confusion Matrix (rows=true, cols=predicted):")
    classes = ['U', 'M', 'C', 'X']
    cm = build_confusion_matrix(true_labels, pred_labels, classes)
    print(cm.to_string())

    # Metrics for UeccDNA
    print("\n" + "=" * 70)
    print("UeccDNA Classification Metrics:")
    u_metrics = compute_metrics(true_labels, pred_labels, 'U')
    print(f"  TP (True U → Pred U): {u_metrics['TP']:,}")
    print(f"  FP (True !U → Pred U): {u_metrics['FP']:,}")
    print(f"  FN (True U → Pred !U): {u_metrics['FN']:,}")
    print(f"  TN (True !U → Pred !U): {u_metrics['TN']:,}")
    print(f"  Precision: {u_metrics['Precision']:.4f}")
    print(f"  Recall: {u_metrics['Recall']:.4f}")
    print(f"  F1 Score: {u_metrics['F1']:.4f}")

    # Metrics for MeccDNA
    print("\n" + "=" * 70)
    print("MeccDNA Classification Metrics:")
    m_metrics = compute_metrics(true_labels, pred_labels, 'M')
    print(f"  TP (True M → Pred M): {m_metrics['TP']:,}")
    print(f"  FP (True !M → Pred M): {m_metrics['FP']:,}")
    print(f"  FN (True M → Pred !M): {m_metrics['FN']:,}")
    print(f"  TN (True !M → Pred !M): {m_metrics['TN']:,}")
    print(f"  Precision: {m_metrics['Precision']:.4f}")
    print(f"  Recall: {m_metrics['Recall']:.4f}")
    print(f"  F1 Score: {m_metrics['F1']:.4f}")

    # Error analysis
    print("\n" + "=" * 70)
    print("Error Analysis:")
    errors = analyze_errors(read_to_type, read_to_pred)
    print(f"  U → M (UeccDNA misclassified as MeccDNA): {len(errors['U_as_M']):,}")
    print(f"  M → U (MeccDNA misclassified as UeccDNA): {len(errors['M_as_U']):,}")
    print(f"  U → X (UeccDNA unclassified): {len(errors['U_as_X']):,}")
    print(f"  M → X (MeccDNA unclassified): {len(errors['M_as_X']):,}")
    print(f"  C → U (CeccDNA misclassified as UeccDNA): {len(errors['C_as_U']):,}")
    print(f"  C → M (CeccDNA misclassified as MeccDNA): {len(errors['C_as_M']):,}")

    # Sample error cases
    print("\n" + "=" * 70)
    print("Sample Error Cases (first 5 each):")
    for error_type, reads in errors.items():
        if reads:
            print(f"\n  {error_type}:")
            for r in reads[:5]:
                print(f"    - {r}")


if __name__ == '__main__':
    main()
