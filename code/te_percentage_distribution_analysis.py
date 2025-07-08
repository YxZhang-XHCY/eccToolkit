#!/usr/bin/env python3
"""
Single & Compound Transposon Percentage Analysis Script
Analyze motif_percent and anno_Percent distributions across different samples and classes
Using 5 bins: 0-20%, 20-40%, 40-60%, 60-80%, 80-100% (including >100%)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time

# Remove any font settings related to Chinese
# rcParams['font.sans-serif'] = ['SimHei', ...]
# rcParams['axes.unicode_minus'] = False


def create_percentage_bins(values):
    """
    Create percentage bins: 0-20%, 20-40%, 40-60%, 60-80%, 80-100% (including >100%)
    """
    bins = [0, 20, 40, 60, 80, float("inf")]
    labels = ["0-20%", "20-40%", "40-60%", "60-80%", "80-100%"]
    values_capped = np.where(values > 100, 100, values)
    return pd.cut(
        values_capped, bins=bins, labels=labels, include_lowest=True, right=False
    )


def analyze_single_compound_te(input_file):
    """
    Main analysis for single/compound transposon percentages
    """
    start_time = time.time()
    print("Reading input data...")
    df = pd.read_csv(input_file)
    print(f"Number of records: {len(df):,}")

    print("\n=== Basic Stats ===")
    print(f"Number of samples: {df['Sample'].nunique()}")
    print(f"Number of TE classes: {df['Class'].nunique()}")
    print(f"Sample list: {', '.join(sorted(df['Sample'].unique()))}")
    print(f"TE class list: {', '.join(sorted(df['Class'].unique()))}")

    print("\n=== Data Quality ===")
    motif_na = df["motif_percent"].isna().sum()
    anno_na = df["anno_Percent"].isna().sum()
    print(f"motif_percent NA: {motif_na:,} ({motif_na/len(df)*100:.1f}%)")
    print(f"anno_Percent NA: {anno_na:,} ({anno_na/len(df)*100:.1f}%)")

    df_motif = df.dropna(subset=["motif_percent"]).copy()
    df_anno = df.dropna(subset=["anno_Percent"]).copy()
    print(f"motif_percent used: {len(df_motif):,}")
    print(f"anno_Percent used: {len(df_anno):,}")

    print("\n=== Binning Percentages ===")
    df_motif["motif_bin"] = create_percentage_bins(df_motif["motif_percent"])
    df_anno["anno_bin"] = create_percentage_bins(df_anno["anno_Percent"])

    motif_over100 = (df_motif["motif_percent"] > 100).sum()
    anno_over100 = (df_anno["anno_Percent"] > 100).sum()
    print(f"motif_percent > 100%: {motif_over100:,}")
    print(f"anno_Percent > 100%: {anno_over100:,}")

    print("\n" + "=" * 50)
    print("MOTIF_PERCENT Distribution Analysis")
    print("=" * 50)
    analyze_percentage_distribution(df_motif, "motif_percent", "motif_bin", "motif")

    print("\n" + "=" * 50)
    print("ANNO_PERCENT Distribution Analysis")
    print("=" * 50)
    analyze_percentage_distribution(df_anno, "anno_Percent", "anno_bin", "anno")

    create_visualizations(df_motif, df_anno)
    generate_detailed_tables(df_motif, df_anno)

    end_time = time.time()
    print(f"\n✅ Analysis complete! Time elapsed: {end_time - start_time:.2f} seconds")


def analyze_percentage_distribution(df, percent_col, bin_col, analysis_type):
    """
    Analyze percentage distributions
    """
    print(f"\n--- {percent_col.upper()} Stats ---")
    stats = df[percent_col].describe()
    print(f"Mean: {stats['mean']:.2f}%")
    print(f"Median: {stats['50%']:.2f}%")
    print(f"Std: {stats['std']:.2f}%")
    print(f"Min: {stats['min']:.2f}%")
    print(f"Max: {stats['max']:.2f}%")

    print(f"\n--- {percent_col.upper()} Bin Counts ---")
    bin_counts = df[bin_col].value_counts().sort_index()
    bin_percent = (bin_counts / len(df) * 100).round(1)
    for bin_name, count in bin_counts.items():
        pct = bin_percent[bin_name]
        print(f"{bin_name}: {count:,} records ({pct}%)")

    print(f"\n--- {percent_col.upper()} by Sample ---")
    sample_stats = df.groupby(["Sample", bin_col]).size().unstack(fill_value=0)
    sample_percent = sample_stats.div(sample_stats.sum(axis=1), axis=0) * 100
    print("Sample bin distribution (percent):")
    print(sample_percent.round(1))

    class_stats = df.groupby(["Class", bin_col]).size().unstack(fill_value=0)
    class_percent = class_stats.div(class_stats.sum(axis=1), axis=0) * 100
    print(f"\nClass bin distribution (percent):")
    print(class_percent.round(1))

    print(f"\n--- {percent_col.upper()} by Sample x Class ---")
    combined_stats = (
        df.groupby(["Sample", "Class", bin_col]).size().unstack(fill_value=0)
    )
    combined_percent = combined_stats.div(combined_stats.sum(axis=1), axis=0) * 100
    print("Sample x Class bin distribution (percent):")
    print(combined_percent.round(1))

    sample_percent.to_csv(f"{analysis_type}_percent_by_sample.csv")
    class_percent.to_csv(f"{analysis_type}_percent_by_class.csv")
    combined_percent.to_csv(f"{analysis_type}_percent_by_sample_class.csv")

    print(f"\n📊 Stats tables saved:")
    print(f"  - {analysis_type}_percent_by_sample.csv")
    print(f"  - {analysis_type}_percent_by_class.csv")
    print(f"  - {analysis_type}_percent_by_sample_class.csv")


def create_visualizations(df_motif, df_anno):
    """
    Generate figures (all English, no Chinese)
    """
    print("\n=== Creating plots ===")
    plt.style.use("default")
    sns.set_palette("husl")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(
        "Distribution of Single/Compound Transposon Percentage",
        fontsize=16,
        fontweight="bold",
    )

    axes[0, 0].hist(
        df_motif["motif_percent"],
        bins=50,
        alpha=0.7,
        color="skyblue",
        edgecolor="black",
    )
    axes[0, 0].set_title("Motif Percent Histogram")
    axes[0, 0].set_xlabel("Motif Percent (%)")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].axvline(
        df_motif["motif_percent"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {df_motif["motif_percent"].mean():.1f}%',
    )
    axes[0, 0].legend()

    motif_bin_counts = df_motif["motif_bin"].value_counts().sort_index()
    axes[0, 1].bar(
        range(len(motif_bin_counts)),
        motif_bin_counts.values,
        color="lightgreen",
        edgecolor="black",
    )
    axes[0, 1].set_title("Motif Percent Bin Distribution")
    axes[0, 1].set_xlabel("Bin")
    axes[0, 1].set_ylabel("Count")
    axes[0, 1].set_xticks(range(len(motif_bin_counts)))
    axes[0, 1].set_xticklabels(motif_bin_counts.index, rotation=45)

    sample_motif = (
        df_motif.groupby(["Sample", "motif_bin"]).size().unstack(fill_value=0)
    )
    sample_motif_pct = sample_motif.div(sample_motif.sum(axis=1), axis=0) * 100
    sample_motif_pct.plot(kind="bar", stacked=True, ax=axes[0, 2], colormap="viridis")
    axes[0, 2].set_title("Motif Percent by Sample")
    axes[0, 2].set_xlabel("Sample")
    axes[0, 2].set_ylabel("Percent (%)")
    axes[0, 2].legend(title="Bin", bbox_to_anchor=(1.05, 1), loc="upper left")
    axes[0, 2].tick_params(axis="x", rotation=45)

    axes[1, 0].hist(
        df_anno["anno_Percent"], bins=50, alpha=0.7, color="orange", edgecolor="black"
    )
    axes[1, 0].set_title("Anno Percent Histogram")
    axes[1, 0].set_xlabel("Anno Percent (%)")
    axes[1, 0].set_ylabel("Frequency")
    axes[1, 0].axvline(
        df_anno["anno_Percent"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {df_anno["anno_Percent"].mean():.1f}%',
    )
    axes[1, 0].legend()

    anno_bin_counts = df_anno["anno_bin"].value_counts().sort_index()
    axes[1, 1].bar(
        range(len(anno_bin_counts)),
        anno_bin_counts.values,
        color="salmon",
        edgecolor="black",
    )
    axes[1, 1].set_title("Anno Percent Bin Distribution")
    axes[1, 1].set_xlabel("Bin")
    axes[1, 1].set_ylabel("Count")
    axes[1, 1].set_xticks(range(len(anno_bin_counts)))
    axes[1, 1].set_xticklabels(anno_bin_counts.index, rotation=45)

    sample_anno = df_anno.groupby(["Sample", "anno_bin"]).size().unstack(fill_value=0)
    sample_anno_pct = sample_anno.div(sample_anno.sum(axis=1), axis=0) * 100
    sample_anno_pct.plot(kind="bar", stacked=True, ax=axes[1, 2], colormap="plasma")
    axes[1, 2].set_title("Anno Percent by Sample")
    axes[1, 2].set_xlabel("Sample")
    axes[1, 2].set_ylabel("Percent (%)")
    axes[1, 2].legend(title="Bin", bbox_to_anchor=(1.05, 1), loc="upper left")
    axes[1, 2].tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig("single_compound_te_analysis.png", dpi=300, bbox_inches="tight")
    plt.savefig("single_compound_te_analysis.pdf", bbox_inches="tight")
    print("📈 Plots saved: single_compound_te_analysis.png/pdf")
    plt.show()


def generate_detailed_tables(df_motif, df_anno):
    """
    Generate summary tables (English headers)
    """
    print("\n=== Generating summary tables ===")
    summary_stats = []

    summary_stats.append(
        {
            "Type": "motif_percent",
            "Sample": "ALL",
            "TE_Class": "ALL",
            "Count": len(df_motif),
            "Mean": df_motif["motif_percent"].mean(),
            "Median": df_motif["motif_percent"].median(),
            "Std": df_motif["motif_percent"].std(),
            "0-20%": (df_motif["motif_bin"] == "0-20%").sum(),
            "20-40%": (df_motif["motif_bin"] == "20-40%").sum(),
            "40-60%": (df_motif["motif_bin"] == "40-60%").sum(),
            "60-80%": (df_motif["motif_bin"] == "60-80%").sum(),
            "80-100%": (df_motif["motif_bin"] == "80-100%").sum(),
        }
    )
    summary_stats.append(
        {
            "Type": "anno_Percent",
            "Sample": "ALL",
            "TE_Class": "ALL",
            "Count": len(df_anno),
            "Mean": df_anno["anno_Percent"].mean(),
            "Median": df_anno["anno_Percent"].median(),
            "Std": df_anno["anno_Percent"].std(),
            "0-20%": (df_anno["anno_bin"] == "0-20%").sum(),
            "20-40%": (df_anno["anno_bin"] == "20-40%").sum(),
            "40-60%": (df_anno["anno_bin"] == "40-60%").sum(),
            "60-80%": (df_anno["anno_bin"] == "60-80%").sum(),
            "80-100%": (df_anno["anno_bin"] == "80-100%").sum(),
        }
    )
    for sample in df_motif["Sample"].unique():
        sample_motif = df_motif[df_motif["Sample"] == sample]
        summary_stats.append(
            {
                "Type": "motif_percent",
                "Sample": sample,
                "TE_Class": "ALL",
                "Count": len(sample_motif),
                "Mean": sample_motif["motif_percent"].mean(),
                "Median": sample_motif["motif_percent"].median(),
                "Std": sample_motif["motif_percent"].std(),
                "0-20%": (sample_motif["motif_bin"] == "0-20%").sum(),
                "20-40%": (sample_motif["motif_bin"] == "20-40%").sum(),
                "40-60%": (sample_motif["motif_bin"] == "40-60%").sum(),
                "60-80%": (sample_motif["motif_bin"] == "60-80%").sum(),
                "80-100%": (sample_motif["motif_bin"] == "80-100%").sum(),
            }
        )
        sample_anno = df_anno[df_anno["Sample"] == sample]
        if len(sample_anno) > 0:
            summary_stats.append(
                {
                    "Type": "anno_Percent",
                    "Sample": sample,
                    "TE_Class": "ALL",
                    "Count": len(sample_anno),
                    "Mean": sample_anno["anno_Percent"].mean(),
                    "Median": sample_anno["anno_Percent"].median(),
                    "Std": sample_anno["anno_Percent"].std(),
                    "0-20%": (sample_anno["anno_bin"] == "0-20%").sum(),
                    "20-40%": (sample_anno["anno_bin"] == "20-40%").sum(),
                    "40-60%": (sample_anno["anno_bin"] == "40-60%").sum(),
                    "60-80%": (sample_anno["anno_bin"] == "60-80%").sum(),
                    "80-100%": (sample_anno["anno_bin"] == "80-100%").sum(),
                }
            )
    summary_df = pd.DataFrame(summary_stats)
    summary_df = summary_df.round(2)
    summary_df.to_csv("single_compound_te_summary.csv", index=False)
    print("📋 Summary table saved: single_compound_te_summary.csv")
    print("\nPreview:")
    print(summary_df.head(10))


def main():
    input_file = "TE.single_compound.csv"
    try:
        analyze_single_compound_te(input_file)
        print(f"\n✅ Analysis completed!")
        print("📁 Output files:")
        print("  - motif_percent_by_sample.csv")
        print("  - motif_percent_by_class.csv")
        print("  - motif_percent_by_sample_class.csv")
        print("  - anno_percent_by_sample.csv")
        print("  - anno_percent_by_class.csv")
        print("  - anno_percent_by_sample_class.csv")
        print("  - single_compound_te_summary.csv")
        print("  - single_compound_te_analysis.png/pdf")
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}")
        print("Please make sure TE.single_compound.csv exists in the current directory")
    except Exception as e:
        print(f"❌ An error occurred: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
