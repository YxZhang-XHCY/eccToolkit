#!/usr/bin/env python3
"""
eccDNA-SHARP: eccDNA Hotspot Analysis with Refined Precision
精确边界的eccDNA热点检测分析流程
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import argparse
import os
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class eccDNAHotspotAnalyzer:
    """eccDNA热点精确边界检测分析器"""
    
    def __init__(self):
        self.window_data = {}
        self.raw_eccdna = {}
        self.chromosomes = []
        self.samples = []
        
    def load_window_data(self, file_10kb, file_50kb, file_100kb):
        """加载多尺度窗口数据"""
        print("加载多尺度窗口数据...")
        self.window_data['10kb'] = pd.read_csv(file_10kb)
        self.window_data['50kb'] = pd.read_csv(file_50kb)
        self.window_data['100kb'] = pd.read_csv(file_100kb)
        
        # 获取样本名称
        self.samples = [col for col in self.window_data['10kb'].columns 
                       if col not in ['Chrom', 'WindowStart', 'WindowEnd']]
        
        # 获取染色体列表
        self.chromosomes = sorted(self.window_data['10kb']['Chrom'].unique())
        
        print(f"检测到 {len(self.samples)} 个样本: {', '.join(self.samples)}")
        print(f"检测到 {len(self.chromosomes)} 条染色体")
        
    def load_raw_eccdna(self, raw_files):
        """加载原始eccDNA数据"""
        print("\n加载原始eccDNA数据...")
        for i, file in enumerate(raw_files):
            sample_name = self.samples[i] if i < len(self.samples) else f"sample_{i}"
            df = pd.read_csv(file)
            self.raw_eccdna[sample_name] = df
            print(f"  {sample_name}: {len(df)} 个eccDNA")
            
    def identify_candidate_hotspots(self, chrom, threshold_factor=2.0):
        """使用多尺度窗口分析识别候选热点区域"""
        candidates = []
        
        for scale, df in self.window_data.items():
            chrom_data = df[df['Chrom'] == chrom].copy()
            
            # 计算每个窗口的平均信号和变异系数
            signal_cols = [col for col in chrom_data.columns 
                          if col not in ['Chrom', 'WindowStart', 'WindowEnd']]
            
            chrom_data['mean_signal'] = chrom_data[signal_cols].mean(axis=1)
            chrom_data['cv'] = chrom_data[signal_cols].std(axis=1) / (chrom_data['mean_signal'] + 1e-6)
            
            # 计算背景水平和阈值
            background = np.median(chrom_data['mean_signal'])
            mad = np.median(np.abs(chrom_data['mean_signal'] - background))
            threshold = background + threshold_factor * mad * 1.4826  # MAD to std
            
            # 识别热点窗口
            hotspot_windows = chrom_data[chrom_data['mean_signal'] > threshold]
            
            # 合并相邻窗口
            if len(hotspot_windows) > 0:
                current_start = None
                current_end = None
                
                for _, window in hotspot_windows.iterrows():
                    if current_start is None:
                        current_start = window['WindowStart']
                        current_end = window['WindowEnd']
                    elif window['WindowStart'] == current_end:
                        current_end = window['WindowEnd']
                    else:
                        candidates.append({
                            'chrom': chrom,
                            'start': current_start,
                            'end': current_end,
                            'scale': scale,
                            'mean_signal': chrom_data[
                                (chrom_data['WindowStart'] >= current_start) & 
                                (chrom_data['WindowEnd'] <= current_end)
                            ]['mean_signal'].mean()
                        })
                        current_start = window['WindowStart']
                        current_end = window['WindowEnd']
                
                if current_start is not None:
                    candidates.append({
                        'chrom': chrom,
                        'start': current_start,
                        'end': current_end,
                        'scale': scale,
                        'mean_signal': chrom_data[
                            (chrom_data['WindowStart'] >= current_start) & 
                            (chrom_data['WindowEnd'] <= current_end)
                        ]['mean_signal'].mean()
                    })
        
        return candidates
    
    def refine_hotspot_boundaries(self, chrom, start, end, resolution=1000):
        """使用原始eccDNA数据精确定义热点边界"""
        refined_regions = []
        
        # 扩展搜索区域
        search_start = max(0, start - 50000)
        search_end = end + 50000
        
        # 创建高分辨率信号轨道
        bins = np.arange(search_start, search_end + resolution, resolution)
        signal_tracks = {}
        
        for sample, df in self.raw_eccdna.items():
            # 筛选在搜索区域内的eccDNA
            region_eccdna = df[
                (df['eChr'] == chrom) & 
                (df['eEnd'] >= search_start) & 
                (df['eStart'] <= search_end)
            ]
            
            # 计算每个bin的信号强度
            signal_array = np.zeros(len(bins) - 1)
            for _, ecc in region_eccdna.iterrows():
                # 找到eccDNA覆盖的bins
                start_bin = np.searchsorted(bins, ecc['eStart'])
                end_bin = np.searchsorted(bins, ecc['eEnd'])
                if start_bin < len(signal_array) and end_bin > 0:
                    signal_array[max(0, start_bin):min(len(signal_array), end_bin)] += 1
            
            signal_tracks[sample] = signal_array
        
        # 计算平均信号
        mean_signal = np.mean(list(signal_tracks.values()), axis=0)
        
        # 平滑信号
        smoothed_signal = gaussian_filter1d(mean_signal, sigma=3)
        
        # 检测信号边界
        boundaries = self.detect_signal_boundaries(smoothed_signal, bins[:-1])
        
        # 验证边界区域
        for b_start, b_end in boundaries:
            if b_start <= end and b_end >= start:  # 与原始候选区域重叠
                # 计算区域统计
                region_stats = self.calculate_region_statistics(
                    chrom, b_start, b_end, signal_tracks, bins[:-1]
                )
                
                if region_stats['mean_signal'] > 0 and region_stats['consistency'] > 0.5:
                    refined_regions.append({
                        'chrom': chrom,
                        'start': b_start,
                        'end': b_end,
                        'length': b_end - b_start,
                        **region_stats
                    })
        
        return refined_regions
    
    def detect_signal_boundaries(self, signal_array, positions, min_peak_height=0.5):
        """检测信号边界"""
        # 注意这里专门引入 signal 模块，避免变量冲突
        from scipy import signal as sp_signal
        boundaries = []
        
        # 计算背景噪声水平
        background = np.median(signal_array)
        noise_level = np.median(np.abs(signal_array - background)) * 1.4826
        threshold = background + 2 * noise_level
        
        # 找到信号峰
        peaks, properties = sp_signal.find_peaks(signal_array, height=max(threshold, min_peak_height))
        
        if len(peaks) == 0:
            return boundaries
        
        # 对每个峰找到边界
        for peak in peaks:
            # 向左找边界
            left_boundary = peak
            while left_boundary > 0 and signal_array[left_boundary] > threshold * 0.5:
                left_boundary -= 1
            
            # 向右找边界
            right_boundary = peak
            while right_boundary < len(signal_array) - 1 and signal_array[right_boundary] > threshold * 0.5:
                right_boundary += 1
            
            # 转换为基因组坐标
            start_pos = int(positions[left_boundary])
            end_pos = int(positions[min(right_boundary + 1, len(positions) - 1)])
            
            # 合并重叠区域
            if boundaries and start_pos <= boundaries[-1][1]:
                boundaries[-1] = (boundaries[-1][0], max(end_pos, boundaries[-1][1]))
            else:
                boundaries.append((start_pos, end_pos))
        
        return boundaries
    
    def calculate_region_statistics(self, chrom, start, end, signal_tracks, positions):
        """计算区域统计信息"""
        # 找到区域对应的bins
        start_idx = np.searchsorted(positions, start)
        end_idx = np.searchsorted(positions, end)
        
        # 计算每个样本的信号
        sample_signals = []
        for sample, signal_array in signal_tracks.items():
            region_signal = signal_array[start_idx:end_idx]
            sample_signals.append(np.sum(region_signal))
        
        # 计算统计信息
        stats = {
            'mean_signal': np.mean(sample_signals),
            'max_signal': np.max(sample_signals),
            'min_signal': np.min(sample_signals),
            'consistency': np.sum(np.array(sample_signals) > 0) / len(sample_signals),
            'cv': np.std(sample_signals) / (np.mean(sample_signals) + 1e-6)
        }
        
        # 计算富集倍数
        total_eccdna = sum(len(df) for df in self.raw_eccdna.values())
        region_length = end - start
        genome_length = 3e9  # 假设基因组大小
        expected = total_eccdna * region_length / genome_length
        stats['enrichment'] = stats['mean_signal'] / (expected / len(self.samples) + 1e-6)
        
        return stats
    
    def classify_hotspots(self, hotspots):
        """对热点进行分类"""
        for hotspot in hotspots:
            consistency = hotspot['consistency']
            enrichment = hotspot['enrichment']
            cv = hotspot['cv']
            
            if consistency >= 0.8 and enrichment >= 5:
                hotspot['category'] = 'core'  # 核心热点
            elif consistency >= 0.6 and enrichment >= 3:
                hotspot['category'] = 'recurrent'  # 复现热点
            else:
                hotspot['category'] = 'variable'  # 可变热点
        
        return hotspots
    
    def permutation_test_boundaries(self, hotspot, n_permutations=1000):
        """对热点边界进行置换检验"""
        chrom = hotspot['chrom']
        start = hotspot['start']
        end = hotspot['end']
        observed_score = hotspot['mean_signal']
        
        # 收集该染色体上的所有eccDNA
        all_eccdna = []
        for sample, df in self.raw_eccdna.items():
            chrom_eccdna = df[df['eChr'] == chrom]
            all_eccdna.extend(chrom_eccdna[['eStart', 'eEnd']].values.tolist())
        
        # 置换检验
        permuted_scores = []
        chrom_length = self.window_data['10kb'][
            self.window_data['10kb']['Chrom'] == chrom
        ]['WindowEnd'].max()
        
        for _ in range(n_permutations):
            # 随机移动热点位置
            shift = np.random.randint(-chrom_length//2, chrom_length//2)
            perm_start = max(0, start + shift)
            perm_end = min(chrom_length, end + shift)
            
            # 计算置换区域的信号
            perm_signal = 0
            for ecc_start, ecc_end in all_eccdna:
                if ecc_end >= perm_start and ecc_start <= perm_end:
                    perm_signal += 1
            
            permuted_scores.append(perm_signal / len(self.samples))
        
        # 计算p值
        p_value = np.sum(np.array(permuted_scores) >= observed_score) / n_permutations
        hotspot['boundary_pvalue'] = p_value
        
        return hotspot
    
    def analyze(self, output_prefix, do_permutation=True):
        """执行完整的分析流程"""
        all_hotspots = []
        
        print("\n开始热点检测分析...")
        for chrom in self.chromosomes:
            print(f"\n处理染色体 {chrom}...")
            
            # 步骤1：识别候选热点
            candidates = self.identify_candidate_hotspots(chrom)
            print(f"  识别到 {len(candidates)} 个候选区域")
            
            # 步骤2：精确定义边界
            for candidate in candidates:
                refined = self.refine_hotspot_boundaries(
                    chrom, 
                    candidate['start'], 
                    candidate['end']
                )
                
                for hotspot in refined:
                    # 步骤3：置换检验（可选）
                    if do_permutation:
                        hotspot = self.permutation_test_boundaries(hotspot)
                    
                    all_hotspots.append(hotspot)
            
            print(f"  精确定义了 {len([h for h in all_hotspots if h['chrom'] == chrom])} 个热点")
        
        # 步骤4：热点分类
        all_hotspots = self.classify_hotspots(all_hotspots)
        
        # 步骤5：输出结果
        self.save_results(all_hotspots, output_prefix)
        
        return all_hotspots
    
    def save_results(self, hotspots, output_prefix):
        """保存分析结果"""
        # 创建输出目录
        output_dir = os.path.dirname(output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 保存热点列表
        df_hotspots = pd.DataFrame(hotspots)
        df_hotspots = df_hotspots.sort_values(['chrom', 'start'])
        df_hotspots.to_csv(f"{output_prefix}_hotspots.csv", index=False)
        
        # 保存BED格式（用于基因组浏览器）
        with open(f"{output_prefix}_hotspots.bed", 'w') as f:
            for _, hotspot in df_hotspots.iterrows():
                score = min(1000, int(hotspot['enrichment'] * 100))
                f.write(f"{hotspot['chrom']}\t{hotspot['start']}\t{hotspot['end']}\t"
                       f"hotspot_{hotspot['category']}\t{score}\t.\n")
        
        # 保存统计摘要
        with open(f"{output_prefix}_summary.txt", 'w') as f:
            f.write("eccDNA热点精确边界检测分析结果摘要\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"总热点数量: {len(hotspots)}\n")
            f.write(f"核心热点: {len([h for h in hotspots if h['category'] == 'core'])}\n")
            f.write(f"复现热点: {len([h for h in hotspots if h['category'] == 'recurrent'])}\n")
            f.write(f"可变热点: {len([h for h in hotspots if h['category'] == 'variable'])}\n\n")
            
            if len(hotspots) > 0:
                f.write(f"平均热点长度: {np.mean([h['length'] for h in hotspots]):.1f} bp\n")
                f.write(f"最小热点长度: {min([h['length'] for h in hotspots])} bp\n")
                f.write(f"最大热点长度: {max([h['length'] for h in hotspots])} bp\n")
                f.write(f"平均富集倍数: {np.mean([h['enrichment'] for h in hotspots]):.2f}\n")
                
                if 'boundary_pvalue' in hotspots[0]:
                    sig_hotspots = [h for h in hotspots if h.get('boundary_pvalue', 1) < 0.05]
                    f.write(f"\n边界显著的热点数量 (p<0.05): {len(sig_hotspots)}\n")
        
        print(f"\n结果已保存至:")
        print(f"  - {output_prefix}_hotspots.csv")
        print(f"  - {output_prefix}_hotspots.bed")
        print(f"  - {output_prefix}_summary.txt")


def main():
    parser = argparse.ArgumentParser(
        description='eccDNA-SHARP: eccDNA热点精确边界检测分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python eccDNA-SHARP.py \\
    --window-10kb data/U87_window_10kb_matrix.csv \\
    --window-50kb data/U87_window_50kb_matrix.csv \\
    --window-100kb data/U87_window_100kb_matrix.csv \\
    --raw-files data/U87_rep*.csv \\
    --output-prefix results/U87
        """
    )
    
    parser.add_argument('--window-10kb', required=True, 
                       help='10kb窗口统计文件')
    parser.add_argument('--window-50kb', required=True, 
                       help='50kb窗口统计文件')
    parser.add_argument('--window-100kb', required=True, 
                       help='100kb窗口统计文件')
    parser.add_argument('--raw-files', nargs='+', required=True,
                       help='原始eccDNA文件列表（CSV格式）')
    parser.add_argument('--output-prefix', required=True,
                       help='输出文件前缀')
    parser.add_argument('--threshold-factor', type=float, default=2.0,
                       help='热点检测阈值因子（默认: 2.0）')
    parser.add_argument('--resolution', type=int, default=1000,
                       help='精细边界检测分辨率（bp，默认: 1000）')
    parser.add_argument('--min-hotspot-size', type=int, default=1000,
                       help='最小热点大小（bp，默认: 1000）')
    parser.add_argument('--permutation-test', action='store_true',
                       help='对热点边界进行置换检验')
    parser.add_argument('--n-permutations', type=int, default=1000,
                       help='置换检验次数（默认: 1000）')
    
    args = parser.parse_args()
    
    # 创建分析器
    analyzer = eccDNAHotspotAnalyzer()
    
    # 加载数据
    analyzer.load_window_data(
        args.window_10kb,
        args.window_50kb,
        args.window_100kb
    )
    analyzer.load_raw_eccdna(args.raw_files)
    
    # 执行分析
    hotspots = analyzer.analyze(
        args.output_prefix,
        do_permutation=args.permutation_test
    )
    
    print(f"\n分析完成！共检测到 {len(hotspots)} 个精确边界的热点区域。")


if __name__ == "__main__":
    main()
