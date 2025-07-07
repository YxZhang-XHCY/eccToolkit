#!/usr/bin/env python3
"""
Analysis of anno_Percent Distribution by Class (MeccDNA vs UeccDNA)
- Only analyze anno_Percent
- Separate analysis for each class
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time

def create_percentage_bins(values):
    bins = [0, 20, 40, 60, 80, float('inf')]
    labels = ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
    values_capped = np.where(values > 100, 100, values)
    return pd.cut(values_capped, bins=bins, labels=labels, include_lowest=True, right=False)

def analyze_anno_percent_by_class(df, te_class):
    print(f"\n=== Analysis for Class: {te_class} ===")
    df_class = df[df['Class'] == te_class].dropna(subset=['anno_Percent']).copy()
    print(f"Number of records: {len(df_class)}")

    # Binning
    df_class['anno_bin'] = create_percentage_bins(df_class['anno_Percent'])
    over100 = (df_class['anno_Percent'] > 100).sum()
    print(f"anno_Percent > 100%: {over100}")

    # Stats
    stats = df_class['anno_Percent'].describe()
    print(f"Mean: {stats['mean']:.2f}%")
    print(f"Median: {stats['50%']:.2f}%")
    print(f"Std: {stats['std']:.2f}%")
    print(f"Min: {stats['min']:.2f}%")
    print(f"Max: {stats['max']:.2f}%")

    # Bin distribution
    bin_counts = df_class['anno_bin'].value_counts().sort_index()
    bin_percent = (bin_counts / len(df_class) * 100).round(1)
    print("\nBin distribution:")
    for bin_name, count in bin_counts.items():
        pct = bin_percent[bin_name]
        print(f"{bin_name}: {count} ({pct}%)")

    # By Sample
    sample_stats = df_class.groupby(['Sample', 'anno_bin'], observed=True).size().unstack(fill_value=0)
    sample_percent = sample_stats.div(sample_stats.sum(axis=1), axis=0) * 100

    # Output tables
    sample_percent.to_csv(f'anno_percent_by_sample_{te_class}.csv')
    bin_counts.to_csv(f'anno_bin_count_{te_class}.csv', header=['count'])
    print(f"\nTables saved: anno_percent_by_sample_{te_class}.csv, anno_bin_count_{te_class}.csv")

    # Visualization
    plot_anno_percent(df_class, sample_percent, te_class)

    # Summary table
    summary = {
        'Class': te_class,
        'Count': len(df_class),
        'Mean': stats['mean'],
        'Median': stats['50%'],
        'Std': stats['std'],
        'Min': stats['min'],
        'Max': stats['max'],
        '0-20%': (df_class['anno_bin'] == '0-20%').sum(),
        '20-40%': (df_class['anno_bin'] == '20-40%').sum(),
        '40-60%': (df_class['anno_bin'] == '40-60%').sum(),
        '60-80%': (df_class['anno_bin'] == '60-80%').sum(),
        '80-100%': (df_class['anno_bin'] == '80-100%').sum()
    }
    return summary

def plot_anno_percent(df_class, sample_percent, te_class):
    plt.style.use('default')
    sns.set_palette("husl")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'anno_Percent Distribution - {te_class}', fontsize=15, fontweight='bold')

    # Histogram
    axes[0].hist(df_class['anno_Percent'], bins=50, alpha=0.7, color='orange', edgecolor='black')
    axes[0].set_title('anno_Percent Histogram')
    axes[0].set_xlabel('anno_Percent (%)')
    axes[0].set_ylabel('Frequency')
    axes[0].axvline(df_class['anno_Percent'].mean(), color='red', linestyle='--', label=f'Mean: {df_class["anno_Percent"].mean():.1f}%')
    axes[0].legend()

    # Bin bar
    bin_counts = df_class['anno_bin'].value_counts().sort_index()
    axes[1].bar(range(len(bin_counts)), bin_counts.values, color='salmon', edgecolor='black')
    axes[1].set_title('anno_Percent Bin Count')
    axes[1].set_xlabel('Bin')
    axes[1].set_ylabel('Count')
    axes[1].set_xticks(range(len(bin_counts)))
    axes[1].set_xticklabels(bin_counts.index, rotation=45)

    # Stacked bar by sample
    sample_percent.plot(kind='bar', stacked=True, ax=axes[2], colormap='plasma')
    axes[2].set_title('anno_Percent by Sample (Stacked)')
    axes[2].set_xlabel('Sample')
    axes[2].set_ylabel('Percent (%)')
    axes[2].legend(title='Bin', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(f'anno_percent_{te_class}_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'anno_percent_{te_class}_analysis.pdf', bbox_inches='tight')
    print(f"Plots saved: anno_percent_{te_class}_analysis.png/pdf")
    plt.close(fig)

def main():
    input_file = "TE.single_compound.csv"
    try:
        start_time = time.time()
        print("Reading input data...")
        df = pd.read_csv(input_file)
        print(f"Total records: {len(df)}")
        class_list = df['Class'].dropna().unique()
        summary_all = []
        for te_class in class_list:
            summary = analyze_anno_percent_by_class(df, te_class)
            summary_all.append(summary)
        # Save summary
        summary_df = pd.DataFrame(summary_all)
        summary_df = summary_df.round(2)
        summary_df.to_csv('anno_percent_by_class_summary.csv', index=False)
        print("\nSummary saved: anno_percent_by_class_summary.csv")
        print(summary_df)
        print(f"\n✅ Analysis completed in {time.time()-start_time:.2f}s!")
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}")
        print("Please make sure TE.single_compound.csv exists in the current directory")
    except Exception as e:
        print(f"❌ An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
