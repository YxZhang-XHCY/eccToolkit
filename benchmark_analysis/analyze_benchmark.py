#!/usr/bin/env python3
"""Benchmark analysis: compare detection results with ground truth."""

import pandas as pd
from collections import defaultdict

def parse_region(region_str):
    """Parse region string like 'Chr1:100-200' to (chrom, start, end)."""
    if ':' not in region_str or '-' not in region_str:
        return None
    # Handle format like 'Chr1:100-200'
    parts = region_str.rsplit(':', 1)  # Split from right to handle complex chrom names
    if len(parts) != 2:
        return None
    chrom, pos = parts
    pos_parts = pos.split('-')
    if len(pos_parts) != 2:
        return None
    try:
        start, end = int(pos_parts[0]), int(pos_parts[1])
        return (chrom, start, end)
    except ValueError:
        return None

def parse_regions(regions_str):
    """Parse multiple regions separated by ';' or '|'."""
    regions = []
    # Handle both ';' and '|' separators
    import re
    parts = re.split(r'[;|]', regions_str)
    for r in parts:
        parsed = parse_region(r.strip())
        if parsed:
            regions.append(parsed)
    return regions

def overlap_ratio(r1, r2):
    """Calculate overlap ratio between two regions."""
    chrom1, s1, e1 = r1
    chrom2, s2, e2 = r2
    if chrom1 != chrom2:
        return 0.0
    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    if overlap_start >= overlap_end:
        return 0.0
    overlap_len = overlap_end - overlap_start
    len1 = e1 - s1
    len2 = e2 - s2
    return overlap_len / min(len1, len2)

def load_truth_bed(filepath):
    """Load ground truth BED file."""
    truth = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            fragments = parts[7] if len(parts) > 7 else None
            truth[name] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'fragments': fragments,
                'matched': False
            }
    return truth

def load_detection_results(filepath):
    """Load detection results CSV."""
    df = pd.read_csv(filepath)
    results = []
    for _, row in df.iterrows():
        regions = parse_regions(row['Regions'])
        results.append({
            'id': row['eccDNA_id'],
            'regions': regions,
            'type': row['eccDNA_type'],
            'state': row['State'],
            'length': row['Length'],
            'matched': False,
            'matched_truth': None
        })
    return results

def match_simple(truth_dict, detections, min_overlap=0.5, multi_region=False):
    """Match simple eccDNA (UeccDNA/MeccDNA) based on region overlap.

    For MeccDNA (multi_region=True), any of the detection regions can match the truth.
    """
    tp = 0
    matched_truth_ids = set()

    for det in detections:
        det_regions = det['regions']
        if not multi_region and len(det_regions) != 1:
            continue

        best_match = None
        best_overlap = 0

        for truth_id, truth_info in truth_dict.items():
            if truth_info['matched']:
                continue
            truth_region = (truth_info['chrom'], truth_info['start'], truth_info['end'])

            # For multi-region, check all regions
            for det_region in det_regions:
                ov = overlap_ratio(det_region, truth_region)
                if ov > best_overlap and ov >= min_overlap:
                    best_overlap = ov
                    best_match = truth_id

        if best_match:
            tp += 1
            det['matched'] = True
            det['matched_truth'] = best_match
            truth_dict[best_match]['matched'] = True
            matched_truth_ids.add(best_match)

    return tp, matched_truth_ids

def match_chimeric(truth_dict, detections, min_overlap=0.3):
    """Match chimeric eccDNA based on fragment overlap."""
    tp = 0
    matched_truth_ids = set()

    for det in detections:
        det_regions = det['regions']
        if len(det_regions) < 2:
            continue

        best_match = None
        best_score = 0

        for truth_id, truth_info in truth_dict.items():
            if truth_info['matched']:
                continue
            if not truth_info['fragments']:
                continue

            # Parse truth fragments
            truth_frags = parse_regions(truth_info['fragments'])
            if len(truth_frags) < 2:
                continue

            # Calculate match score
            matched_frags = 0
            for det_r in det_regions:
                for truth_r in truth_frags:
                    if overlap_ratio(det_r, truth_r) >= min_overlap:
                        matched_frags += 1
                        break

            score = matched_frags / max(len(det_regions), len(truth_frags))
            if score > best_score and score >= 0.5:
                best_score = score
                best_match = truth_id

        if best_match:
            tp += 1
            det['matched'] = True
            det['matched_truth'] = best_match
            truth_dict[best_match]['matched'] = True
            matched_truth_ids.add(best_match)

    return tp, matched_truth_ids

def main():
    # Load data
    print("=" * 60)
    print("eccDNA Detection Benchmark Analysis")
    print("=" * 60)

    truth_unique = load_truth_bed('Demo_1400.unique.bed')
    truth_multi = load_truth_bed('Demo_1400.multi.bed')
    truth_chimeric = load_truth_bed('Demo_1400.chimeric.bed')

    detections = load_detection_results('sample_merged_output.csv')

    # Separate detections by type
    det_unique = [d for d in detections if d['type'] == 'UeccDNA']
    det_multi = [d for d in detections if d['type'] == 'MeccDNA']
    det_chimeric = [d for d in detections if d['type'] == 'CeccDNA']

    print(f"\nGround Truth:")
    print(f"  UeccDNA: {len(truth_unique)}")
    print(f"  MeccDNA: {len(truth_multi)}")
    print(f"  CeccDNA: {len(truth_chimeric)}")

    print(f"\nDetection Results:")
    print(f"  UeccDNA: {len(det_unique)}")
    print(f"  MeccDNA: {len(det_multi)}")
    print(f"  CeccDNA: {len(det_chimeric)}")

    # Match UeccDNA
    tp_u, _ = match_simple(truth_unique, det_unique, min_overlap=0.5, multi_region=False)
    fn_u = len(truth_unique) - tp_u
    fp_u = len(det_unique) - tp_u

    # Match MeccDNA (multi-region: any of the detection regions can match)
    tp_m, _ = match_simple(truth_multi, det_multi, min_overlap=0.5, multi_region=True)
    fn_m = len(truth_multi) - tp_m
    fp_m = len(det_multi) - tp_m

    # Match CeccDNA
    tp_c, _ = match_chimeric(truth_chimeric, det_chimeric, min_overlap=0.3)
    fn_c = len(truth_chimeric) - tp_c
    fp_c = len(det_chimeric) - tp_c

    # Calculate metrics
    def calc_metrics(tp, fp, fn):
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        return precision, recall, f1

    print("\n" + "=" * 60)
    print("Benchmark Results")
    print("=" * 60)

    print(f"\n{'Type':<12} {'TP':>6} {'FP':>6} {'FN':>6} {'Precision':>10} {'Recall':>10} {'F1':>10}")
    print("-" * 62)

    p, r, f1 = calc_metrics(tp_u, fp_u, fn_u)
    print(f"{'UeccDNA':<12} {tp_u:>6} {fp_u:>6} {fn_u:>6} {p:>10.4f} {r:>10.4f} {f1:>10.4f}")

    p, r, f1 = calc_metrics(tp_m, fp_m, fn_m)
    print(f"{'MeccDNA':<12} {tp_m:>6} {fp_m:>6} {fn_m:>6} {p:>10.4f} {r:>10.4f} {f1:>10.4f}")

    p, r, f1 = calc_metrics(tp_c, fp_c, fn_c)
    print(f"{'CeccDNA':<12} {tp_c:>6} {fp_c:>6} {fn_c:>6} {p:>10.4f} {r:>10.4f} {f1:>10.4f}")

    # Overall
    tp_all = tp_u + tp_m + tp_c
    fp_all = fp_u + fp_m + fp_c
    fn_all = fn_u + fn_m + fn_c
    p, r, f1 = calc_metrics(tp_all, fp_all, fn_all)
    print("-" * 62)
    print(f"{'Overall':<12} {tp_all:>6} {fp_all:>6} {fn_all:>6} {p:>10.4f} {r:>10.4f} {f1:>10.4f}")

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total Ground Truth: {len(truth_unique) + len(truth_multi) + len(truth_chimeric)}")
    print(f"Total Detections:   {len(detections)}")
    print(f"True Positives:     {tp_all}")
    print(f"False Positives:    {fp_all}")
    print(f"False Negatives:    {fn_all}")

if __name__ == '__main__':
    main()
