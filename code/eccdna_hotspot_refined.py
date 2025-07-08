#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os
from scipy.ndimage import label
from tqdm import tqdm


def load_window_candidates(window_files, threshold_factor=2.0):
    """
    输入多尺度窗口文件，合并所有候选区段。
    返回：候选区段表 DataFrame（Chrom,Start,End）
    """
    candidates = []
    for scale, file in window_files.items():
        df = pd.read_csv(file)
        sample_cols = [
            c for c in df.columns if c not in ["Chrom", "WindowStart", "WindowEnd"]
        ]
        df["mean_signal"] = df[sample_cols].mean(axis=1)
        background = np.median(df["mean_signal"])
        mad = np.median(np.abs(df["mean_signal"] - background))
        threshold = background + threshold_factor * mad * 1.4826
        df = df[df["mean_signal"] > threshold]
        cdf = df[["Chrom", "WindowStart", "WindowEnd"]].copy()
        cdf["scale"] = scale
        candidates.append(cdf)
    all_cand = pd.concat(candidates)
    # 合并重叠区段
    merged = []
    for chrom, group in all_cand.groupby("Chrom"):
        intervals = (
            group[["WindowStart", "WindowEnd"]].sort_values("WindowStart").values
        )
        s, e = None, None
        for interval in intervals:
            if s is None:
                s, e = interval
            elif interval[0] <= e:
                e = max(e, interval[1])
            else:
                merged.append([chrom, s, e])
                s, e = interval
        if s is not None:
            merged.append([chrom, s, e])
    return pd.DataFrame(merged, columns=["Chrom", "Start", "End"])


def load_all_eccdna_beds(raw_files):
    """
    读入所有eccDNA bed/csv，返回每个样本DataFrame列表
    """
    dfs = []
    for f in raw_files:
        df = pd.read_csv(f)
        if not set(["eChr", "eStart", "eEnd"]).issubset(df.columns):
            raise Exception(f"{f} 缺少eChr/eStart/eEnd列")
        dfs.append(df[["eChr", "eStart", "eEnd"]].copy())
    return dfs


def get_highres_signal(chrom, start, end, dfs, bin_size=1000):
    """
    合并所有样本，统计±50kb 1kb bin信号
    返回：bin信号数组、bin起止坐标
    """
    search_start = max(0, start - 50000)
    search_end = end + 50000
    bins = np.arange(search_start, search_end + bin_size, bin_size)
    signal_track = np.zeros(len(bins) - 1)
    for df in dfs:
        df1 = df[
            (df["eChr"] == chrom)
            & (df["eEnd"] >= search_start)
            & (df["eStart"] <= search_end)
        ]
        for _, row in df1.iterrows():
            sidx = np.searchsorted(bins, row["eStart"], side="right") - 1
            eidx = np.searchsorted(bins, row["eEnd"], side="right") - 1
            sidx = max(0, sidx)
            eidx = min(len(signal_track) - 1, eidx)
            signal_track[sidx : eidx + 1] += 1
    return signal_track, bins


def find_hotspot_regions(signal_track, bins, orig_start, orig_end, min_size=1000):
    """
    用信号阈值法和最大连续区段法识别热点区
    返回：区间(start,end)列表（仅与候选区overlap的部分）
    """
    background = np.median(signal_track)
    noise = np.median(np.abs(signal_track - background)) * 1.4826
    threshold = background + 2 * noise
    is_hot = signal_track > threshold
    lbl, num = label(is_hot)
    hotspots = []
    for i in range(1, num + 1):
        idxs = np.where(lbl == i)[0]
        region_start = bins[idxs[0]]
        region_end = bins[idxs[-1] + 1]
        if region_end - region_start >= min_size:
            # 必须和原粗定位区段有overlap
            if region_end > orig_start and region_start < orig_end:
                hotspots.append((region_start, region_end))
    return hotspots


def compute_hotspot_stats(chrom, s, e, dfs, genome_size=3e9):
    """
    信号强度、重现性、富集倍数
    """
    signals = []
    for df in dfs:
        df1 = df[(df["eChr"] == chrom) & (df["eEnd"] >= s) & (df["eStart"] <= e)]
        signals.append(len(df1))
    mean_signal = np.mean(signals)
    consistency = np.sum(np.array(signals) > 0) / len(signals)
    total_ecc = sum([len(df) for df in dfs])
    region_len = e - s
    expected = total_ecc * region_len / genome_size
    enrichment = mean_signal / (expected / len(dfs) + 1e-9)
    return mean_signal, enrichment, consistency, min(signals), max(signals)


def permutation_test(chrom, s, e, dfs, chrom_length, n_perm=1000):
    """
    置换热点区间并统计信号，返回p值
    """
    obs = np.mean(
        [
            len(df[(df["eChr"] == chrom) & (df["eEnd"] >= s) & (df["eStart"] <= e)])
            for df in dfs
        ]
    )
    region_len = e - s
    null = []
    for _ in range(n_perm):
        shift = np.random.randint(0, chrom_length - region_len)
        ps = shift
        pe = shift + region_len
        perm_val = np.mean(
            [
                len(
                    df[
                        (df["eChr"] == chrom)
                        & (df["eEnd"] >= ps)
                        & (df["eStart"] <= pe)
                    ]
                )
                for df in dfs
            ]
        )
        null.append(perm_val)
    null = np.array(null)
    pval = (np.sum(null >= obs) + 1) / (n_perm + 1)
    return pval


def classify_hotspot(consistency, enrichment):
    """
    热点分类
    """
    if consistency >= 0.8 and enrichment >= 5:
        return "core"
    elif consistency >= 0.6 and enrichment >= 3:
        return "recurrent"
    else:
        return "variable"


def main():
    parser = argparse.ArgumentParser(description="eccDNA热点精确边界检测")
    parser.add_argument("--window-10kb", required=True)
    parser.add_argument("--window-50kb", required=True)
    parser.add_argument("--window-100kb", required=True)
    parser.add_argument("--raw-files", nargs="+", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--threshold-factor", type=float, default=2.0)
    parser.add_argument("--min-hotspot-size", type=int, default=1000)
    parser.add_argument("--permutation-test", action="store_true")
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--bin-size", type=int, default=1000)
    args = parser.parse_args()

    window_files = {
        "10kb": args.window_10kb,
        "50kb": args.window_50kb,
        "100kb": args.window_100kb,
    }
    # Step1: 粗定位
    print("候选热点区段粗定位...")
    candidates = load_window_candidates(window_files, args.threshold_factor)
    print(f"候选热点区段数: {len(candidates)}")

    # Step2: 高分辨率精细化
    print("加载eccDNA原始文件...")
    dfs = load_all_eccdna_beds(args.raw_files)
    results = []
    print("热点区段精确边界检测...")
    for idx, row in tqdm(candidates.iterrows(), total=len(candidates)):
        chrom, start, end = row["Chrom"], int(row["Start"]), int(row["End"])
        # 读取染色体最大长度（这里用窗口文件推断）
        chrom_length = 0
        for win in window_files.values():
            dfc = pd.read_csv(win)
            dfc_chrom = dfc[dfc["Chrom"] == chrom]
            chrom_length = max(chrom_length, dfc_chrom["WindowEnd"].max())
        signal_track, bins = get_highres_signal(chrom, start, end, dfs, args.bin_size)
        regions = find_hotspot_regions(
            signal_track, bins, start, end, args.min_hotspot_size
        )
        for s, e in regions:
            mean_signal, enrichment, consistency, min_sig, max_sig = (
                compute_hotspot_stats(chrom, s, e, dfs)
            )
            pval = None
            if args.permutation_test:
                pval = permutation_test(
                    chrom, s, e, dfs, chrom_length, args.n_permutations
                )
            results.append(
                {
                    "chrom": chrom,
                    "start": int(s),
                    "end": int(e),
                    "length": int(e - s),
                    "mean_signal": mean_signal,
                    "enrichment": enrichment,
                    "consistency": consistency,
                    "min_signal": min_sig,
                    "max_signal": max_sig,
                    "category": classify_hotspot(consistency, enrichment),
                    "boundary_pvalue": pval,
                }
            )
    # 输出结果
    dfres = pd.DataFrame(results)
    if not os.path.exists(os.path.dirname(args.output_prefix)) and os.path.dirname(
        args.output_prefix
    ):
        os.makedirs(os.path.dirname(args.output_prefix))
    dfres.to_csv(f"{args.output_prefix}_hotspots.csv", index=False)
    # BED输出
    with open(f"{args.output_prefix}_hotspots.bed", "w") as f:
        for _, row in dfres.iterrows():
            score = min(1000, int(row["enrichment"] * 100))
            f.write(
                f'{row["chrom"]}\t{int(row["start"])}\t{int(row["end"])}\thotspot_{row["category"]}\t{score}\t.\n'
            )
    print("分析完成。结果见：")
    print(f"  {args.output_prefix}_hotspots.csv")
    print(f"  {args.output_prefix}_hotspots.bed")


if __name__ == "__main__":
    main()
