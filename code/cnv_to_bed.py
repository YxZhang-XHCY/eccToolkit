#!/usr/bin/env python3
"""
Split a 1-based CNV.tsv into per-cell-line BED files.

BED output columns:
  1  chrom
  2  start (0-based, inclusive)
  3  end   (1-based, exclusive—unchanged)
  4  SegmentMean
  5  NumProbes
  6  Status
  7  ModelID
  8  ProfileID
"""

import argparse, os
import pandas as pd

def main():
    p = argparse.ArgumentParser(
        description="Convert CNV.tsv to per-cell-line BED files")
    p.add_argument("-i", "--input",  required=True, help="CNV.tsv")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # load tab-separated CNV file (header present)
    df = pd.read_csv(args.input, sep="\t")

    # convert 1-based Start → 0-based for BED
    df["Start0"] = df["Start"] - 1

    bed_cols = ["Chromosome", "Start0", "End",
                "SegmentMean", "NumProbes", "Status",
                "ModelID", "ProfileID"]

    # write one BED per Cell_Lines
    for cell, sub in df.groupby("Cell_Lines", sort=False):
        out_path = os.path.join(args.outdir, f"{cell}.cnv.bed")
        sub[bed_cols].to_csv(out_path,
                             sep="\t",
                             header=False,
                             index=False)

if __name__ == "__main__":
    main()
