#!/usr/bin/env python3
"""
eccDNA富集/耗竭与RNA-seq差异表达基因关联分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom, fisher_exact, chi2_contingency, spearmanr
from matplotlib_venn import venn2, venn3
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings("ignore")

# 设置中文字体
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def load_eccdna_data(filepath):
    """加载eccDNA富集数据"""
    df = pd.read_csv(filepath)
    return df


def load_degs_data(filepath):
    """加载DEGs数据"""
    # 根据您的数据格式，这是一个制表符分隔的文件，没有列名
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["gene", "ensembl_id", "tumor_exp", "normal_exp", "log2fc", "qvalue"],
    )
    return df


def analyze_eccdna_degs_overlap(eccdna_df, degs_df, sample_name=None, fc_threshold=1):
    """分析单个样本的eccDNA与DEGs重叠"""

    # 如果指定了样本名，过滤该样本的数据
    if sample_name:
        eccdna_sample = eccdna_df[eccdna_df["Sample"] == sample_name].copy()
    else:
        eccdna_sample = eccdna_df.copy()

    # 获取显著的eccDNA基因
    eccdna_enriched = set(
        eccdna_sample[
            (eccdna_sample["direction"] == "enrichment")
            & (eccdna_sample["significant"] == True)
        ]["gene"]
    )
    eccdna_depleted = set(
        eccdna_sample[
            (eccdna_sample["direction"] == "depletion")
            & (eccdna_sample["significant"] == True)
        ]["gene"]
    )

    # 获取DEGs（使用指定的FC阈值）
    degs_up = set(degs_df[degs_df["log2fc"] > fc_threshold]["gene"])
    degs_down = set(degs_df[degs_df["log2fc"] < -fc_threshold]["gene"])

    # 计算交集
    results = {
        "enriched_up": eccdna_enriched & degs_up,
        "enriched_down": eccdna_enriched & degs_down,
        "depleted_up": eccdna_depleted & degs_up,
        "depleted_down": eccdna_depleted & degs_down,
    }

    # 统计信息
    stats = {
        "n_eccdna_enriched": len(eccdna_enriched),
        "n_eccdna_depleted": len(eccdna_depleted),
        "n_degs_up": len(degs_up),
        "n_degs_down": len(degs_down),
        "n_enriched_up": len(results["enriched_up"]),
        "n_enriched_down": len(results["enriched_down"]),
        "n_depleted_up": len(results["depleted_up"]),
        "n_depleted_down": len(results["depleted_down"]),
    }

    return results, stats, eccdna_enriched, eccdna_depleted, degs_up, degs_down


def hypergeometric_test(overlap, set1_size, set2_size, total_genes):
    """超几何检验"""
    # 计算p值
    p_value = hypergeom.sf(overlap - 1, total_genes, set1_size, set2_size)
    return p_value


def fisher_test(overlap, set1_size, set2_size, total_genes):
    """Fisher精确检验"""
    # 构建2x2列联表
    table = [
        [overlap, set1_size - overlap],
        [set2_size - overlap, total_genes - set1_size - set2_size + overlap],
    ]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")
    return odds_ratio, p_value


def perform_gradient_fc_analysis(eccdna_df, degs_df, fc_thresholds=[1, 2, 4, 6]):
    """执行梯度FC阈值分析"""
    results_summary = []

    # 计算背景基因数
    all_genes = set(eccdna_df["gene"]) | set(degs_df["gene"])
    total_genes = len(all_genes)

    # 获取eccDNA基因集（这些不会随FC阈值改变）
    eccdna_enriched_all = set(
        eccdna_df[
            (eccdna_df["direction"] == "enrichment")
            & (eccdna_df["significant"] == True)
        ]["gene"]
    )
    eccdna_depleted_all = set(
        eccdna_df[
            (eccdna_df["direction"] == "depletion") & (eccdna_df["significant"] == True)
        ]["gene"]
    )

    print("\n=== Gradient FC Threshold Analysis ===")
    print(f"FC thresholds: {fc_thresholds}")
    print(f"Total background genes: {total_genes}")
    print(f"eccDNA enriched genes: {len(eccdna_enriched_all)}")
    print(f"eccDNA depleted genes: {len(eccdna_depleted_all)}")

    for fc in fc_thresholds:
        print(f"\n--- FC threshold: {fc} ---")

        # 分析当前FC阈值
        results, stats, _, _, _, _ = analyze_eccdna_degs_overlap(
            eccdna_df, degs_df, fc_threshold=fc
        )

        # 计算统计显著性
        analysis_results = {
            "fc_threshold": fc,
            "n_degs_up": stats["n_degs_up"],
            "n_degs_down": stats["n_degs_down"],
            "n_enriched_up": stats["n_enriched_up"],
            "n_enriched_down": stats["n_enriched_down"],
            "n_depleted_up": stats["n_depleted_up"],
            "n_depleted_down": stats["n_depleted_down"],
        }

        # Fisher检验：enriched vs up
        if stats["n_eccdna_enriched"] > 0 and stats["n_degs_up"] > 0:
            or_eu, p_eu = fisher_test(
                stats["n_enriched_up"],
                stats["n_eccdna_enriched"],
                stats["n_degs_up"],
                total_genes,
            )
            analysis_results["enriched_up_OR"] = or_eu
            analysis_results["enriched_up_pvalue"] = p_eu
            print(
                f"  Enriched ∩ Up: n={stats['n_enriched_up']}, OR={or_eu:.3f}, p={p_eu:.3e}"
            )

        # Fisher检验：depleted vs down
        if stats["n_eccdna_depleted"] > 0 and stats["n_degs_down"] > 0:
            or_dd, p_dd = fisher_test(
                stats["n_depleted_down"],
                stats["n_eccdna_depleted"],
                stats["n_degs_down"],
                total_genes,
            )
            analysis_results["depleted_down_OR"] = or_dd
            analysis_results["depleted_down_pvalue"] = p_dd
            print(
                f"  Depleted ∩ Down: n={stats['n_depleted_down']}, OR={or_dd:.3f}, p={p_dd:.3e}"
            )

        # Fisher检验：enriched vs down (反向)
        if stats["n_eccdna_enriched"] > 0 and stats["n_degs_down"] > 0:
            or_ed, p_ed = fisher_test(
                stats["n_enriched_down"],
                stats["n_eccdna_enriched"],
                stats["n_degs_down"],
                total_genes,
            )
            analysis_results["enriched_down_OR"] = or_ed
            analysis_results["enriched_down_pvalue"] = p_ed
            print(
                f"  Enriched ∩ Down: n={stats['n_enriched_down']}, OR={or_ed:.3f}, p={p_ed:.3e}"
            )

        # Fisher检验：depleted vs up (反向)
        if stats["n_eccdna_depleted"] > 0 and stats["n_degs_up"] > 0:
            or_du, p_du = fisher_test(
                stats["n_depleted_up"],
                stats["n_eccdna_depleted"],
                stats["n_degs_up"],
                total_genes,
            )
            analysis_results["depleted_up_OR"] = or_du
            analysis_results["depleted_up_pvalue"] = p_du
            print(
                f"  Depleted ∩ Up: n={stats['n_depleted_up']}, OR={or_du:.3f}, p={p_du:.3e}"
            )

        results_summary.append(analysis_results)

    return pd.DataFrame(results_summary)


def plot_gradient_analysis(gradient_df):
    """绘制梯度分析结果"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("eccDNA-DEG Association Across FC Thresholds", fontsize=16)

    # 1. 重叠基因数量随FC阈值变化
    ax = axes[0, 0]
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_enriched_up"],
        "o-",
        label="Enriched ∩ Up",
        color="red",
    )
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_depleted_down"],
        "s-",
        label="Depleted ∩ Down",
        color="blue",
    )
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_enriched_down"],
        "^-",
        label="Enriched ∩ Down",
        color="orange",
    )
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_depleted_up"],
        "v-",
        label="Depleted ∩ Up",
        color="green",
    )
    ax.set_xlabel("FC Threshold")
    ax.set_ylabel("Number of Overlapping Genes")
    ax.set_title("Overlapping Genes vs FC Threshold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Odds Ratio随FC阈值变化
    ax = axes[0, 1]
    if "enriched_up_OR" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            gradient_df["enriched_up_OR"],
            "o-",
            label="Enriched-Up",
            color="red",
        )
    if "depleted_down_OR" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            gradient_df["depleted_down_OR"],
            "s-",
            label="Depleted-Down",
            color="blue",
        )
    if "enriched_down_OR" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            gradient_df["enriched_down_OR"],
            "^-",
            label="Enriched-Down",
            color="orange",
        )
    if "depleted_up_OR" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            gradient_df["depleted_up_OR"],
            "v-",
            label="Depleted-Up",
            color="green",
        )
    ax.axhline(y=1, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("FC Threshold")
    ax.set_ylabel("Odds Ratio")
    ax.set_title("Association Strength (OR) vs FC Threshold")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale("log")

    # 3. P值随FC阈值变化
    ax = axes[1, 0]
    if "enriched_up_pvalue" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            -np.log10(gradient_df["enriched_up_pvalue"]),
            "o-",
            label="Enriched-Up",
            color="red",
        )
    if "depleted_down_pvalue" in gradient_df.columns:
        ax.plot(
            gradient_df["fc_threshold"],
            -np.log10(gradient_df["depleted_down_pvalue"]),
            "s-",
            label="Depleted-Down",
            color="blue",
        )
    ax.axhline(
        y=-np.log10(0.05), color="black", linestyle="--", alpha=0.5, label="p=0.05"
    )
    ax.set_xlabel("FC Threshold")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Statistical Significance vs FC Threshold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. DEG数量随FC阈值变化
    ax = axes[1, 1]
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_degs_up"],
        "o-",
        label="Upregulated DEGs",
        color="red",
    )
    ax.plot(
        gradient_df["fc_threshold"],
        gradient_df["n_degs_down"],
        "s-",
        label="Downregulated DEGs",
        color="blue",
    )
    ax.set_xlabel("FC Threshold")
    ax.set_ylabel("Number of DEGs")
    ax.set_title("DEG Count vs FC Threshold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_fc_specific_venn(eccdna_df, degs_df, fc_threshold=2):
    """为特定FC阈值绘制韦恩图"""
    results, stats, eccdna_enriched, eccdna_depleted, degs_up, degs_down = (
        analyze_eccdna_degs_overlap(eccdna_df, degs_df, fc_threshold=fc_threshold)
    )

    fig = plot_venn_diagrams(
        eccdna_enriched,
        eccdna_depleted,
        degs_up,
        degs_down,
        f"(FC threshold: {fc_threshold})",
    )
    return fig


def plot_venn_diagrams(
    eccdna_enriched, eccdna_depleted, degs_up, degs_down, sample_name=""
):
    """绘制韦恩图"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"eccDNA vs DEGs Overlap Analysis {sample_name}", fontsize=16)

    # 1. eccDNA富集 vs RNA-seq上调
    ax = axes[0, 0]
    venn2(
        [eccdna_enriched, degs_up], set_labels=["eccDNA Enriched", "RNA-seq Up"], ax=ax
    )
    ax.set_title("eccDNA Enriched vs RNA-seq Upregulated")

    # 2. eccDNA耗竭 vs RNA-seq下调
    ax = axes[0, 1]
    venn2(
        [eccdna_depleted, degs_down],
        set_labels=["eccDNA Depleted", "RNA-seq Down"],
        ax=ax,
    )
    ax.set_title("eccDNA Depleted vs RNA-seq Downregulated")

    # 3. eccDNA富集 vs RNA-seq下调（反向）
    ax = axes[1, 0]
    venn2(
        [eccdna_enriched, degs_down],
        set_labels=["eccDNA Enriched", "RNA-seq Down"],
        ax=ax,
    )
    ax.set_title("eccDNA Enriched vs RNA-seq Downregulated (Inverse)")

    # 4. eccDNA耗竭 vs RNA-seq上调（反向）
    ax = axes[1, 1]
    venn2(
        [eccdna_depleted, degs_up], set_labels=["eccDNA Depleted", "RNA-seq Up"], ax=ax
    )
    ax.set_title("eccDNA Depleted vs RNA-seq Upregulated (Inverse)")

    plt.tight_layout()
    return fig


def analyze_sample_overlap(eccdna_df):
    """分析不同样本间的eccDNA模式重叠"""
    samples = eccdna_df["Sample"].unique()

    # 存储每个样本的基因集
    sample_enriched = {}
    sample_depleted = {}

    for sample in samples:
        sample_data = eccdna_df[eccdna_df["Sample"] == sample]
        sample_enriched[sample] = set(
            sample_data[
                (sample_data["direction"] == "enrichment")
                & (sample_data["significant"] == True)
            ]["gene"]
        )
        sample_depleted[sample] = set(
            sample_data[
                (sample_data["direction"] == "depletion")
                & (sample_data["significant"] == True)
            ]["gene"]
        )

    return sample_enriched, sample_depleted


def plot_sample_overlap_venn(sample_genes, title, max_samples=3):
    """绘制样本间的韦恩图"""
    samples = list(sample_genes.keys())[:max_samples]

    if len(samples) == 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        venn2(
            [sample_genes[samples[0]], sample_genes[samples[1]]],
            set_labels=samples,
            ax=ax,
        )
        ax.set_title(title)
    elif len(samples) >= 3:
        fig, ax = plt.subplots(figsize=(8, 6))
        venn3(
            [
                sample_genes[samples[0]],
                sample_genes[samples[1]],
                sample_genes[samples[2]],
            ],
            set_labels=samples[:3],
            ax=ax,
        )
        ax.set_title(title)

    return fig


def create_overlap_heatmap(eccdna_df, degs_df):
    """创建基因在不同分析中的状态热图"""
    # 获取所有样本
    samples = eccdna_df["Sample"].unique()

    # 获取所有相关基因
    all_genes = set()
    for sample in samples:
        sample_data = eccdna_df[eccdna_df["Sample"] == sample]
        significant_genes = sample_data[sample_data["significant"] == True]["gene"]
        all_genes.update(significant_genes)

    # 添加DEGs中的基因
    all_genes.update(degs_df["gene"])

    # 创建基因状态矩阵
    gene_list = sorted(list(all_genes))
    status_matrix = pd.DataFrame(index=gene_list)

    # 添加每个样本的eccDNA状态
    for sample in samples:
        sample_data = eccdna_df[eccdna_df["Sample"] == sample]
        enriched = sample_data[
            (sample_data["direction"] == "enrichment")
            & (sample_data["significant"] == True)
        ]["gene"]
        depleted = sample_data[
            (sample_data["direction"] == "depletion")
            & (sample_data["significant"] == True)
        ]["gene"]

        status_matrix[f"{sample}_eccDNA"] = 0
        status_matrix.loc[status_matrix.index.isin(enriched), f"{sample}_eccDNA"] = 1
        status_matrix.loc[status_matrix.index.isin(depleted), f"{sample}_eccDNA"] = -1

    # 添加DEG状态
    degs_dict = dict(zip(degs_df["gene"], degs_df["log2fc"]))
    status_matrix["DEG_status"] = 0
    for gene in gene_list:
        if gene in degs_dict:
            if degs_dict[gene] > 0:
                status_matrix.loc[gene, "DEG_status"] = 1
            else:
                status_matrix.loc[gene, "DEG_status"] = -1

    return status_matrix


def perform_chi_square_trend_test(gradient_results):
    """执行卡方趋势检验，评估关联强度随FC阈值的变化趋势"""
    from scipy.stats import chi2

    print("\n=== Chi-square Trend Test ===")

    # 对每种关联类型进行趋势检验
    associations = ["enriched_up", "depleted_down", "enriched_down", "depleted_up"]
    trend_results = {}

    for assoc in associations:
        n_overlap = gradient_results[f"n_{assoc}"].values
        fc_thresholds = gradient_results["fc_threshold"].values

        # 计算线性趋势的卡方统计量
        n = len(fc_thresholds)
        x = np.arange(n)

        # 加权线性回归
        weights = n_overlap / n_overlap.sum() if n_overlap.sum() > 0 else np.ones(n) / n
        x_mean = np.sum(x * weights)
        y_mean = np.sum(n_overlap * weights)

        num = np.sum(weights * (x - x_mean) * (n_overlap - y_mean))
        den = np.sum(weights * (x - x_mean) ** 2)

        if den > 0:
            slope = num / den
            chi2_stat = slope**2 * den
            p_value = 1 - chi2.cdf(chi2_stat, df=1)

            trend_results[assoc] = {
                "slope": slope,
                "chi2": chi2_stat,
                "p_value": p_value,
                "trend": "increasing" if slope > 0 else "decreasing",
            }

            print(
                f"{assoc}: slope={slope:.3f}, χ²={chi2_stat:.3f}, p={p_value:.3e}, trend={trend_results[assoc]['trend']}"
            )

    return trend_results


def calculate_enrichment_score(eccdna_df, degs_df):
    """计算每个基因的综合富集评分"""
    # 合并eccDNA数据（取平均值）
    eccdna_summary = (
        eccdna_df[eccdna_df["significant"] == True]
        .groupby("gene")
        .agg(
            {
                "fold_change": "mean",
                "z_score": "mean",
                "direction": lambda x: x.mode()[0] if len(x.mode()) > 0 else x.iloc[0],
            }
        )
        .reset_index()
    )

    # 合并DEG数据
    merged = pd.merge(
        eccdna_summary, degs_df[["gene", "log2fc", "qvalue"]], on="gene", how="outer"
    )

    # 计算综合评分
    merged["eccdna_score"] = merged["fold_change"].fillna(1) * np.sign(
        merged["z_score"].fillna(0)
    )
    merged["deg_score"] = merged["log2fc"].fillna(0)
    merged["combined_score"] = merged["eccdna_score"] * merged["deg_score"]

    # 分类基因
    def classify_gene(row):
        if pd.isna(row["fold_change"]) or pd.isna(row["log2fc"]):
            return "Single_method"
        elif row["direction"] == "enrichment" and row["log2fc"] > 0:
            return "Concordant_positive"
        elif row["direction"] == "depletion" and row["log2fc"] < 0:
            return "Concordant_negative"
        else:
            return "Discordant"

    merged["category"] = merged.apply(classify_gene, axis=1)

    return merged


def plot_correlation_analysis(eccdna_df, degs_df):
    """绘制eccDNA fold change与RNA-seq log2FC的相关性分析"""
    # 准备数据
    gene_scores = calculate_enrichment_score(eccdna_df, degs_df)
    gene_scores = gene_scores.dropna(subset=["fold_change", "log2fc"])

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. 散点图：eccDNA FC vs RNA-seq log2FC
    ax = axes[0, 0]
    colors = {"enrichment": "red", "depletion": "blue"}
    for direction, color in colors.items():
        subset = gene_scores[gene_scores["direction"] == direction]
        ax.scatter(
            subset["log2fc"],
            np.log2(subset["fold_change"]),
            c=color,
            alpha=0.6,
            label=f"eccDNA {direction}",
            s=50,
        )

    # 添加相关性
    if len(gene_scores) > 3:
        corr, p_val = spearmanr(
            gene_scores["log2fc"], np.log2(gene_scores["fold_change"])
        )
        ax.text(
            0.05,
            0.95,
            f"Spearman ρ = {corr:.3f}\np = {p_val:.3e}",
            transform=ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("RNA-seq Log2 Fold Change")
    ax.set_ylabel("eccDNA Log2 Fold Change")
    ax.set_title("eccDNA vs RNA-seq Fold Changes")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. 密度图：显示基因分布
    ax = axes[0, 1]
    from scipy.stats import gaussian_kde

    # 计算2D核密度估计
    x = gene_scores["log2fc"].values
    y = np.log2(gene_scores["fold_change"].values)

    # 移除无限值
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) > 10:
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)

        scatter = ax.scatter(x, y, c=z, s=50, cmap="viridis", alpha=0.6)
        plt.colorbar(scatter, ax=ax, label="Density")

    ax.set_xlabel("RNA-seq Log2 Fold Change")
    ax.set_ylabel("eccDNA Log2 Fold Change")
    ax.set_title("Gene Distribution Density")
    ax.grid(True, alpha=0.3)

    # 3. 基因分类统计
    ax = axes[1, 0]
    category_counts = gene_scores["category"].value_counts()
    colors_cat = {
        "Concordant_positive": "darkgreen",
        "Concordant_negative": "darkred",
        "Discordant": "orange",
        "Single_method": "gray",
    }

    bars = ax.bar(
        category_counts.index,
        category_counts.values,
        color=[colors_cat.get(x, "gray") for x in category_counts.index],
    )
    ax.set_xlabel("Gene Category")
    ax.set_ylabel("Number of Genes")
    ax.set_title("Gene Classification by eccDNA-DEG Relationship")
    ax.tick_params(axis="x", rotation=45)

    # 添加数值标签
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{int(height)}",
            ha="center",
            va="bottom",
        )

    # 4. 四象限分析
    ax = axes[1, 1]

    # 定义四象限
    quadrants = {
        "Q1": (gene_scores["log2fc"] > 0) & (gene_scores["eccdna_score"] > 1),
        "Q2": (gene_scores["log2fc"] < 0) & (gene_scores["eccdna_score"] > 1),
        "Q3": (gene_scores["log2fc"] < 0) & (gene_scores["eccdna_score"] < 1),
        "Q4": (gene_scores["log2fc"] > 0) & (gene_scores["eccdna_score"] < 1),
    }

    quadrant_names = [
        "Q1: Up/Enriched",
        "Q2: Down/Enriched",
        "Q3: Down/Depleted",
        "Q4: Up/Depleted",
    ]
    quadrant_counts = [quadrants[q].sum() for q in ["Q1", "Q2", "Q3", "Q4"]]

    ax.pie(
        quadrant_counts,
        labels=quadrant_names,
        autopct="%1.1f%%",
        colors=["green", "orange", "red", "yellow"],
    )
    ax.set_title("Quadrant Analysis of Gene Changes")

    plt.tight_layout()
    return fig, gene_scores


def perform_pathway_enrichment_prep(gene_scores):
    """准备基因列表用于通路富集分析"""
    # 创建不同类别的基因列表
    gene_lists = {
        "concordant_positive": gene_scores[
            gene_scores["category"] == "Concordant_positive"
        ]["gene"].tolist(),
        "concordant_negative": gene_scores[
            gene_scores["category"] == "Concordant_negative"
        ]["gene"].tolist(),
        "discordant": gene_scores[gene_scores["category"] == "Discordant"][
            "gene"
        ].tolist(),
        "eccdna_enriched_only": gene_scores[
            (gene_scores["direction"] == "enrichment")
            & (pd.isna(gene_scores["log2fc"]))
        ]["gene"].tolist(),
        "eccdna_depleted_only": gene_scores[
            (gene_scores["direction"] == "depletion") & (pd.isna(gene_scores["log2fc"]))
        ]["gene"].tolist(),
        "deg_up_only": gene_scores[
            (gene_scores["log2fc"] > 0) & (pd.isna(gene_scores["fold_change"]))
        ]["gene"].tolist(),
        "deg_down_only": gene_scores[
            (gene_scores["log2fc"] < 0) & (pd.isna(gene_scores["fold_change"]))
        ]["gene"].tolist(),
    }

    # 保存基因列表
    with open("gene_lists_for_enrichment.txt", "w") as f:
        for category, genes in gene_lists.items():
            f.write(f"=== {category} ({len(genes)} genes) ===\n")
            for gene in genes:
                f.write(f"{gene}\n")
            f.write("\n")

    return gene_lists


def main():
    # 文件路径
    eccdna_file = "merged_gene_eccdna_enrichment_results.csv"
    degs_file = "DEGs_GBM.GEPIA2.Log2FC_1.qValue_0.01.txt"

    # 加载数据
    print("Loading data...")
    eccdna_df = load_eccdna_data(eccdna_file)
    degs_df = load_degs_data(degs_file)

    # 获取所有样本 - 直接使用'Sample'列名
    samples = eccdna_df["Sample"].unique()
    print(f"Found {len(samples)} samples: {', '.join(samples)}")

    # 1. 执行梯度FC阈值分析
    fc_thresholds = [1, 2, 4, 6]
    gradient_results = perform_gradient_fc_analysis(eccdna_df, degs_df, fc_thresholds)

    # 保存梯度分析结果
    gradient_results.to_csv("gradient_fc_analysis_results.csv", index=False)
    print("\nGradient analysis results saved to: gradient_fc_analysis_results.csv")

    # 绘制梯度分析图
    fig_gradient = plot_gradient_analysis(gradient_results)
    plt.savefig("gradient_fc_analysis.png", dpi=300, bbox_inches="tight")
    print("Gradient analysis plot saved to: gradient_fc_analysis.png")

    # 2. 为关键FC阈值绘制详细韦恩图
    for fc in [1, 2, 4]:
        fig_venn = plot_fc_specific_venn(eccdna_df, degs_df, fc_threshold=fc)
        plt.savefig(f"eccdna_degs_venn_fc{fc}.png", dpi=300, bbox_inches="tight")
        print(f"Venn diagram for FC={fc} saved to: eccdna_degs_venn_fc{fc}.png")

    # 3. 分析整体的eccDNA与DEGs重叠（使用FC=1作为基准）
    print("\n=== Overall Analysis (FC=1) ===")
    results, stats, eccdna_enriched, eccdna_depleted, degs_up, degs_down = (
        analyze_eccdna_degs_overlap(eccdna_df, degs_df, fc_threshold=1)
    )

    # 打印统计结果
    print(f"Total eccDNA enriched genes: {stats['n_eccdna_enriched']}")
    print(f"Total eccDNA depleted genes: {stats['n_eccdna_depleted']}")
    print(f"Total upregulated DEGs: {stats['n_degs_up']}")
    print(f"Total downregulated DEGs: {stats['n_degs_down']}")
    print(f"\nOverlap statistics:")
    print(f"  Enriched ∩ Upregulated: {stats['n_enriched_up']}")
    print(f"  Enriched ∩ Downregulated: {stats['n_enriched_down']}")
    print(f"  Depleted ∩ Upregulated: {stats['n_depleted_up']}")
    print(f"  Depleted ∩ Downregulated: {stats['n_depleted_down']}")

    # 4. 分析样本间的重叠
    print("\n=== Sample Overlap Analysis ===")
    sample_enriched, sample_depleted = analyze_sample_overlap(eccdna_df)

    # 打印样本间重叠统计
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i < j:
                overlap_enrich = len(sample_enriched[s1] & sample_enriched[s2])
                overlap_deplete = len(sample_depleted[s1] & sample_depleted[s2])
                print(f"{s1} vs {s2}:")
                print(f"  Enriched overlap: {overlap_enrich}")
                print(f"  Depleted overlap: {overlap_deplete}")

    # 5. 绘制样本间韦恩图
    if len(samples) >= 2:
        fig2 = plot_sample_overlap_venn(
            sample_enriched, "eccDNA Enriched Genes Overlap", 3
        )
        plt.savefig("eccdna_enriched_sample_venn.png", dpi=300, bbox_inches="tight")

        fig3 = plot_sample_overlap_venn(
            sample_depleted, "eccDNA Depleted Genes Overlap", 3
        )
        plt.savefig("eccdna_depleted_sample_venn.png", dpi=300, bbox_inches="tight")

    # 6. 创建热图显示基因状态
    print("\n=== Creating gene status heatmap ===")
    status_matrix = create_overlap_heatmap(eccdna_df, degs_df)

    # 只显示在至少一个条件下显著的基因
    significant_genes = status_matrix[(status_matrix != 0).any(axis=1)]

    if len(significant_genes) > 0:
        # 限制显示的基因数量
        if len(significant_genes) > 100:
            # 选择变化最大的前100个基因
            gene_variance = significant_genes.var(axis=1)
            top_genes = gene_variance.nlargest(100).index
            significant_genes = significant_genes.loc[top_genes]

        plt.figure(figsize=(12, 10))
        sns.heatmap(
            significant_genes.T,
            cmap="RdBu_r",
            center=0,
            cbar_kws={"label": "Status (-1: Depleted/Down, 0: NS, 1: Enriched/Up)"},
            xticklabels=False,
            yticklabels=True,
        )
        plt.xlabel("Genes")
        plt.title("Gene Status Across eccDNA and RNA-seq Analysis")
        plt.tight_layout()
        plt.savefig("gene_status_heatmap.png", dpi=300, bbox_inches="tight")

    # 7. 输出详细结果文件
    with open("eccdna_degs_detailed_results.txt", "w") as f:
        f.write("=== eccDNA vs DEGs Overlap Analysis Results ===\n\n")

        # 写入梯度分析总结
        f.write("=== Gradient FC Analysis Summary ===\n")
        f.write(gradient_results.to_string())
        f.write("\n\n")

        # 对每个FC阈值输出详细基因列表
        for fc in fc_thresholds:
            f.write(f"\n=== FC Threshold: {fc} ===\n")
            results, stats, _, _, _, _ = analyze_eccdna_degs_overlap(
                eccdna_df, degs_df, fc_threshold=fc
            )

            f.write(
                f"\nEnriched & Upregulated ({len(results['enriched_up'])} genes):\n"
            )
            for gene in sorted(results["enriched_up"])[:20]:  # 显示前20个
                f.write(f"  {gene}\n")

            f.write(
                f"\nDepleted & Downregulated ({len(results['depleted_down'])} genes):\n"
            )
            for gene in sorted(results["depleted_down"])[:20]:
                f.write(f"  {gene}\n")

    # 8. 创建基因FC与eccDNA关系的散点图
    print("\n=== Creating scatter plots for gene-level analysis ===")

    # 合并数据以创建散点图
    eccdna_gene_stats = (
        eccdna_df[eccdna_df["significant"] == True]
        .groupby("gene")
        .agg({"fold_change": "mean", "direction": "first"})
        .reset_index()
    )

    # 合并DEG数据
    merged_data = pd.merge(
        eccdna_gene_stats, degs_df[["gene", "log2fc"]], on="gene", how="inner"
    )

    # 创建散点图
    plt.figure(figsize=(10, 8))
    colors = {"enrichment": "red", "depletion": "blue"}
    for direction, color in colors.items():
        subset = merged_data[merged_data["direction"] == direction]
        plt.scatter(
            subset["log2fc"],
            subset["fold_change"],
            c=color,
            alpha=0.6,
            label=f"eccDNA {direction}",
        )

    plt.axhline(y=1, color="gray", linestyle="--", alpha=0.5)
    plt.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    plt.xlabel("RNA-seq Log2 Fold Change")
    plt.ylabel("eccDNA Fold Change")
    plt.title("Relationship between eccDNA and RNA-seq Changes")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("eccdna_rnaseq_scatter.png", dpi=300, bbox_inches="tight")

    # 9. 进行高级相关性分析
    fig_corr, gene_scores = plot_correlation_analysis(eccdna_df, degs_df)
    plt.savefig("eccdna_degs_correlation_analysis.png", dpi=300, bbox_inches="tight")
    print("Correlation analysis saved to: eccdna_degs_correlation_analysis.png")

    # 10. 执行卡方趋势检验
    trend_results = perform_chi_square_trend_test(gradient_results)

    # 11. 准备通路富集分析的基因列表
    gene_lists = perform_pathway_enrichment_prep(gene_scores)
    print("\nGene lists for pathway enrichment saved to: gene_lists_for_enrichment.txt")

    # 12. 创建综合报告
    with open("comprehensive_analysis_report.txt", "w") as f:
        f.write("=== Comprehensive eccDNA-DEG Analysis Report ===\n\n")

        # 梯度分析总结
        f.write("1. GRADIENT FC ANALYSIS\n")
        f.write("-" * 50 + "\n")
        f.write(gradient_results.to_string())
        f.write("\n\n")

        # 趋势检验结果
        f.write("2. CHI-SQUARE TREND TEST RESULTS\n")
        f.write("-" * 50 + "\n")
        for assoc, results in trend_results.items():
            f.write(
                f"{assoc}: {results['trend']} trend, χ²={results['chi2']:.3f}, p={results['p_value']:.3e}\n"
            )
        f.write("\n")

        # 基因分类统计
        f.write("3. GENE CLASSIFICATION SUMMARY\n")
        f.write("-" * 50 + "\n")
        category_counts = gene_scores["category"].value_counts()
        for cat, count in category_counts.items():
            f.write(f"{cat}: {count} genes ({count/len(gene_scores)*100:.1f}%)\n")
        f.write("\n")

        # 最强关联的基因
        f.write("4. TOP CONCORDANT GENES\n")
        f.write("-" * 50 + "\n")

        # 正向一致的基因
        f.write("\nTop Concordant Positive (eccDNA enriched + RNA upregulated):\n")
        concordant_pos = gene_scores[
            gene_scores["category"] == "Concordant_positive"
        ].nlargest(20, "combined_score")
        for _, row in concordant_pos.iterrows():
            f.write(
                f"  {row['gene']}: eccDNA_FC={row['fold_change']:.2f}, RNA_log2FC={row['log2fc']:.2f}\n"
            )

        # 负向一致的基因
        f.write("\nTop Concordant Negative (eccDNA depleted + RNA downregulated):\n")
        concordant_neg = gene_scores[
            gene_scores["category"] == "Concordant_negative"
        ].nsmallest(20, "combined_score")
        for _, row in concordant_neg.iterrows():
            f.write(
                f"  {row['gene']}: eccDNA_FC={row['fold_change']:.2f}, RNA_log2FC={row['log2fc']:.2f}\n"
            )

    print("\nComprehensive report saved to: comprehensive_analysis_report.txt")
    print("\n=== Analysis Complete! ===")
    print("\nAll output files:")
    print("  - gradient_fc_analysis_results.csv (梯度分析数据)")
    print("  - gradient_fc_analysis.png (梯度分析图)")
    print("  - eccdna_degs_venn_fc*.png (不同FC阈值的韦恩图)")
    print("  - eccdna_*_sample_venn.png (样本间韦恩图)")
    print("  - gene_status_heatmap.png (基因状态热图)")
    print("  - eccdna_rnaseq_scatter.png (散点图)")
    print("  - eccdna_degs_correlation_analysis.png (相关性分析)")
    print("  - gene_lists_for_enrichment.txt (用于通路富集的基因列表)")
    print("  - comprehensive_analysis_report.txt (综合分析报告)")
    print("  - eccdna_degs_detailed_results.txt (详细结果)")


if __name__ == "__main__":
    main()
