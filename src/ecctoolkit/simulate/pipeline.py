"""
统一模拟流水线模块

编排 sim-region（区域生成）和 readsim（读段模拟）的完整模拟流程。
"""

import logging
import os
from pathlib import Path
from typing import List, Optional

from .unified_config import UnifiedSimulateConfig

logger = logging.getLogger(__name__)


class SimulatePipeline:
    """编排完整的 eccDNA 模拟流水线"""

    def __init__(
        self,
        config: UnifiedSimulateConfig,
        output_dir: str,
        reference: str,
        skip_region: bool = False,
        input_eccdna: Optional[str] = None,
        skip_sr: bool = False,
        skip_hifi: bool = False,
        skip_ont: bool = False,
        compress: bool = False,
        verbose: bool = False,
        coverage_list: Optional[List[float]] = None,
    ):
        """
        初始化流水线。

        Args:
            config: 统一配置对象
            output_dir: 输出目录
            reference: 参考基因组 FASTA 路径
            skip_region: 跳过区域生成
            input_eccdna: 已有 eccDNA FASTA 路径（跳过区域生成时使用）
            skip_sr: 跳过短读模拟
            skip_hifi: 跳过 HiFi 模拟
            skip_ont: 跳过 ONT 模拟
            compress: 压缩输出文件
            verbose: 详细日志
            coverage_list: 覆盖度列表（支持多覆盖度模拟）
        """
        self.config = config
        self.output_dir = Path(output_dir)
        self.reference = reference
        self.skip_region = skip_region
        self.input_eccdna = input_eccdna
        self.skip_sr = skip_sr
        self.skip_hifi = skip_hifi
        self.skip_ont = skip_ont
        self.compress = compress
        self.verbose = verbose
        # 覆盖度列表，默认使用配置中的单一覆盖度
        self.coverage_list = coverage_list or [config.readsim.params.meancov]

        # 设置日志
        log_level = logging.DEBUG if verbose else logging.INFO
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
        )

    @staticmethod
    def _fmt_cov(c: float) -> str:
        """格式化覆盖度：整数显示为整数，小数保留一位"""
        return f"{int(c)}" if c == int(c) else f"{c:.1f}"

    def run(self) -> None:
        """执行完整的模拟流水线"""
        logger.info("=" * 60)
        logger.info("eccDNA 模拟流水线开始")
        logger.info("=" * 60)

        # 设置输出目录
        self._setup_output_dirs()

        # 保存使用的配置
        self._save_config()

        # 步骤 1: 区域模拟
        eccdna_fasta = self._run_region_simulation()
        logger.info(f"eccDNA FASTA: {eccdna_fasta}")

        # 步骤 2: 读段模拟（支持多覆盖度）
        if not self.config.readsim.enabled:
            logger.info("跳过读段模拟 (--skip-readsim)")
        else:
            self._run_readsim_multi_coverage(eccdna_fasta)

        logger.info("=" * 60)
        logger.info("模拟流水线完成")
        logger.info(f"输出目录: {self.output_dir}")
        if len(self.coverage_list) > 1:
            cov_dirs = [f"sequencing_{self._fmt_cov(c)}X" for c in self.coverage_list]
            logger.info(f"覆盖度目录: {', '.join(cov_dirs)}")
        logger.info("=" * 60)

    def _setup_output_dirs(self) -> None:
        """创建输出目录结构"""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建子目录
        self.region_dir = self.output_dir / "region"

        if not self.skip_region:
            self.region_dir.mkdir(exist_ok=True)

        # RCA 扩增目录
        if self.config.readsim.enabled:
            self.rca_dir = self.output_dir / "rca"
            self.rca_dir.mkdir(exist_ok=True)

        # 多覆盖度时使用 sequencing_{cov}X 目录，单覆盖度时使用 sequencing
        if self.config.readsim.enabled:
            if len(self.coverage_list) == 1:
                self.sequencing_dir = self.output_dir / "sequencing"
                self.sequencing_dir.mkdir(exist_ok=True)
            else:
                # 多覆盖度模式：为每个覆盖度创建目录
                for cov in self.coverage_list:
                    cov_dir = self.output_dir / f"sequencing_{self._fmt_cov(cov)}X"
                    cov_dir.mkdir(exist_ok=True)
                # sequencing_dir 设为第一个覆盖度目录（用于兼容）
                self.sequencing_dir = self.output_dir / f"sequencing_{self._fmt_cov(self.coverage_list[0])}X"

        self.lite_dir = self.sequencing_dir

    def _save_config(self) -> None:
        """保存使用的配置"""
        config_path = self.output_dir / "config_used.yaml"
        self.config.to_yaml(config_path)
        logger.info(f"配置已保存: {config_path}")

    def _run_region_simulation(self) -> Path:
        """运行区域模拟或返回已有 eccDNA"""
        if self.skip_region:
            if self.input_eccdna:
                logger.info(f"使用已有 eccDNA: {self.input_eccdna}")
                return Path(self.input_eccdna)
            else:
                raise ValueError("跳过区域生成时必须提供 --input-eccdna")

        logger.info("-" * 40)
        logger.info("步骤 1: 区域模拟 (sim-region)")
        logger.info("-" * 40)

        from .region import run_region_simulation

        # 提取配置参数
        region = self.config.region
        ld = region.length_distribution
        cls = region.classification
        mm = region.minimap

        # 直接在 region_dir 下输出，使用 sample 作为文件前缀
        output_prefix = self.config.sample

        run_region_simulation(
            reference=self.reference,
            output_prefix=output_prefix,
            output_dir=str(self.region_dir),
            num_unique=region.num_unique,
            num_multi=region.num_multi,
            num_chimeric=region.num_chimeric,
            threads=self.config.threads,
            seed=self.config.seed or 42,
            # 长度分布
            mode=ld.mode,
            sigma=ld.sigma,
            tail_weight=ld.tail_weight,
            tail_min=ld.tail_min,
            tail_max=ld.tail_max,
            min_length=ld.min_length,
            max_length=ld.max_length,
            max_length_multi=region.max_length_multi,
            # 分类阈值
            identity=cls.identity,
            min_coverage=cls.min_coverage,
            length_consistency=cls.length_consistency,
            multi_coverage=cls.multi_coverage,
            max_secondary=cls.max_secondary,
            no_hit_policy=cls.no_hit_policy,
            # minimap2 设置
            split_by_length=mm.split_by_length,
            split_length=mm.split_length,
            minimap_preset_short=mm.preset_short,
            minimap_preset_long=mm.preset_long,
            # 候选倍数
            candidate_multiplier_u=region.candidate_multiplier_u,
            candidate_multiplier_m=region.candidate_multiplier_m,
            # 输出选项
            keep_tmp=region.keep_tmp,
            verbose=self.verbose,
        )

        return self.region_dir / f"{self.config.sample}.all.fa"

    def _generate_linear_dna(self, pos_csv) -> tuple:
        """
        生成背景线性 DNA，数量与 eccDNA 1:1。

        从参考基因组随机采样区域，长度分布与 eccDNA 相似。

        Returns:
            (bed_rows, csv_rows): 线性 DNA 的 bed 和 csv 数据
        """
        import numpy as np
        from tqdm import tqdm
        from .lite_utils import utilities

        utils = utilities(self.reference)

        # 获取 eccDNA 的长度分布
        ecc_lengths = pos_csv['length'].values
        num_linear = len(ecc_lengths)

        if num_linear == 0:
            return [], []

        # 获取参考基因组的染色体信息
        genome_info = utils.genome_length()
        chroms = genome_info.index.tolist()
        chrom_lengths = genome_info['length'].values

        # 按染色体长度加权采样
        chrom_weights = chrom_lengths / chrom_lengths.sum()

        bed_rows = []
        csv_rows = []

        # 设置随机种子
        if self.config.seed is not None:
            np.random.seed(self.config.seed + 12345)  # 使用不同的种子偏移

        logger.info(f"  生成 {num_linear} 个线性 DNA 区域...")

        for i in tqdm(range(num_linear), desc="生成线性 DNA"):
            # 使用对应 eccDNA 的长度（保持长度分布一致）
            length = int(ecc_lengths[i])

            # 随机选择染色体（加权）
            chrom_idx = np.random.choice(len(chroms), p=chrom_weights)
            chrom = chroms[chrom_idx]

            # 确保长度不超过染色体长度，且至少为 100bp
            max_len = int(chrom_lengths[chrom_idx]) - 1000  # 留一些余量
            if max_len < 100:
                # 染色体太短，尝试其他染色体
                for alt_idx in range(len(chroms)):
                    alt_max = int(chrom_lengths[alt_idx]) - 1000
                    if alt_max >= 100:
                        chrom_idx = alt_idx
                        chrom = chroms[chrom_idx]
                        max_len = alt_max
                        break
                else:
                    # 所有染色体都太短，跳过此线性 DNA
                    logger.debug(f"跳过线性 DNA {i}: 所有染色体都太短")
                    continue

            if length > max_len:
                length = max_len
            if length < 100:
                length = 100

            try:
                # 从参考基因组随机采样
                region = utils.random_region(chrom, length)
                chrom, start, end, frag_len, seq = region
            except (RuntimeError, ValueError) as e:
                # 如果采样失败，尝试其他染色体
                logger.debug(f"采样失败 ({chrom}, {length}): {e}, 尝试其他染色体")
                sampled = False
                for _ in range(10):
                    chrom_idx = np.random.choice(len(chroms), p=chrom_weights)
                    chrom = chroms[chrom_idx]
                    try_len = min(length, max(100, int(chrom_lengths[chrom_idx]) - 1000))
                    if try_len < 100:
                        continue
                    try:
                        region = utils.random_region(chrom, try_len)
                        chrom, start, end, frag_len, seq = region
                        sampled = True
                        break
                    except (RuntimeError, ValueError, IndexError) as e:
                        logger.debug(f"备选染色体 {chrom} 采样失败: {e}")
                        continue
                if not sampled:
                    # 如果还是失败，跳过此线性 DNA
                    logger.debug(f"跳过线性 DNA {i}: 采样失败")
                    continue

            linear_id = f"neg_{i+1:06d}"

            bed_rows.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "length": int(frag_len),
                "id": linear_id,
            })

            csv_rows.append({
                "id": linear_id,
                "fragN": 1,
                "region": f"{chrom}:{start}-{end}",
                "length": int(frag_len),
                "seq": seq,
            })

        return bed_rows, csv_rows

    def _convert_region_format(self, eccdna_fasta: Path) -> None:
        """
        将 sim-region 的输出转换为读段模拟需要的格式。

        sim-region 输出:
            - {sample}.all.fa: FASTA 格式，header 包含区域信息
            - {sample}.all.bed: BED 格式

        读段模拟需要:
            - {sample}.pos.bed: chrom, start, end, length, id
            - {sample}.pos.csv: id, fragN, region, length, seq
            - {sample}.neg.bed/csv: 线性 DNA（这里生成空文件）

        所有 RCA 相关的中间文件都放在 rca_dir 中。
        """
        import gzip
        import pandas as pd
        from Bio import SeqIO

        # RCA 中间文件放在 rca_dir 中
        sample_dir = self.rca_dir
        sample_dir.mkdir(parents=True, exist_ok=True)

        logger.info("转换 sim-region 输出格式...")

        # Prefer sim-region BED (contains chimeric fragments) to avoid losing coordinates.
        bed_candidates = [eccdna_fasta.with_suffix(".bed")]
        if eccdna_fasta.suffix == ".gz":
            bed_candidates.insert(0, eccdna_fasta.with_suffix("").with_suffix(".bed"))
        region_bed_path = next((p for p in bed_candidates if p.exists()), None)

        def _open_fasta(path: Path):
            if str(path).endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "r")

        # Parse FASTA sequences for pos.csv (libsim only needs pos.bed).
        with _open_fasta(eccdna_fasta) as handle:
            seq_by_id = {rec.id: str(rec.seq) for rec in SeqIO.parse(handle, "fasta")}
        if not seq_by_id:
            raise ValueError(f"No records found in eccDNA FASTA: {eccdna_fasta}")

        bed_rows = []
        csv_rows = []

        def _parse_fragment(fragment):
            chrom_part, coords = fragment.split(":", 1)
            start_s, end_s = coords.split("-", 1)
            return chrom_part, int(start_s), int(end_s)

        if region_bed_path is None:
            # Fallback: parse coordinates from FASTA headers.
            # sim-region headers:
            #   >UeccDNA_000001 chr1:100-200 length=... type=U
            #   >CeccDNA_000001 chr1:100-150;chr2:200-260 length=... type=C
            with _open_fasta(eccdna_fasta) as handle:
                records = SeqIO.parse(handle, "fasta")
                for rec in records:
                    parts = rec.description.split()
                    if len(parts) < 2:
                        raise ValueError(
                            f"Cannot infer coordinates from FASTA header (missing region field): {rec.description}"
                        )
                    region_field = parts[1]
                    fragments = [_parse_fragment(s) for s in region_field.split(";") if s and s != "."]
                    if not fragments:
                        raise ValueError(f"Cannot infer fragments from FASTA header: {rec.description}")

                    region_parts = []
                    total_len = 0
                    for chrom, start, end in fragments:
                        frag_len = int(end) - int(start)
                        bed_rows.append(
                            {
                                "chrom": chrom,
                                "start": int(start),
                                "end": int(end),
                                "length": frag_len,
                                "id": rec.id,
                            }
                        )
                        region_parts.append(f"{chrom}:{int(start)}-{int(end)}")
                        total_len += frag_len

                    seq = seq_by_id.get(rec.id, "")
                    csv_rows.append(
                        {
                            "id": rec.id,
                            "fragN": len(fragments),
                            "region": "|".join(region_parts),
                            "length": len(seq) if seq else total_len,
                            "seq": seq,
                        }
                    )
        else:
            # BED columns: chrom start end name length strand type fragments
            region_bed = pd.read_csv(region_bed_path, sep="\t", comment="#", header=None, dtype=str)
            if region_bed.empty:
                raise ValueError(f"sim-region BED is empty: {region_bed_path}")

            if region_bed.shape[1] < 7:
                raise ValueError(
                    f"Unexpected sim-region BED format (need >=7 columns): {region_bed_path}"
                )

            columns = ["chrom", "start", "end", "name", "length", "strand", "type"]
            if region_bed.shape[1] >= 8:
                columns.append("fragments")
            region_bed = region_bed.iloc[:, : len(columns)]
            region_bed.columns = columns

            region_bed["start"] = region_bed["start"].astype(int)
            region_bed["end"] = region_bed["end"].astype(int)
            if "length" in region_bed.columns:
                region_bed["length"] = region_bed["length"].astype(int)

            for _, row in region_bed.iterrows():
                region_id = str(row["name"])
                fragments_field = str(row["fragments"]) if "fragments" in row else "."
                fragments_field = fragments_field.strip()

                fragments = []
                if fragments_field and fragments_field not in {".", "nan", "NaN"}:
                    fragments = [
                        _parse_fragment(s.strip())
                        for s in fragments_field.split(";")
                        if s.strip() and s.strip() != "."
                    ]

                if not fragments:
                    fragments = [(str(row["chrom"]), int(row["start"]), int(row["end"]))]

                region_parts = []
                total_len = 0
                for chrom, start, end in fragments:
                    frag_len = int(end) - int(start)
                    bed_rows.append(
                        {
                            "chrom": chrom,
                            "start": int(start),
                            "end": int(end),
                            "length": frag_len,
                            "id": region_id,
                        }
                    )
                    region_parts.append(f"{chrom}:{int(start)}-{int(end)}")
                    total_len += frag_len

                seq = seq_by_id.get(region_id, "")
                csv_rows.append(
                    {
                        "id": region_id,
                        "fragN": len(fragments),
                        "region": "|".join(region_parts),
                        "length": len(seq) if seq else total_len,
                        "seq": seq,
                    }
                )

        # 创建 DataFrame 并保存
        pos_bed = pd.DataFrame(bed_rows, columns=["chrom", "start", "end", "length", "id"])
        pos_csv = pd.DataFrame(csv_rows, columns=["id", "fragN", "region", "length", "seq"])

        # 保存 pos 文件（环形 DNA）
        pos_bed_path = sample_dir / f"{self.config.sample}.pos.bed"
        pos_csv_path = sample_dir / f"{self.config.sample}.pos.csv"

        pos_bed.to_csv(pos_bed_path, sep="\t", header=False, index=False)
        pos_csv.to_csv(pos_csv_path, sep="\t", index=False)

        # 生成线性 DNA（背景）- 与 eccDNA 数量 1:1
        logger.info("生成背景线性 DNA...")
        neg_bed_rows = []
        neg_csv_rows = []
        reference_path = Path(self.reference) if self.reference else None
        if reference_path is None or not reference_path.exists():
            logger.warning("Reference not found; skipping linear DNA generation.")
        else:
            try:
                neg_bed_rows, neg_csv_rows = self._generate_linear_dna(pos_csv)
            except Exception as exc:
                logger.warning(f"Linear DNA generation failed ({exc}); writing empty neg files.")

        neg_bed_path = sample_dir / f"{self.config.sample}.neg.bed"
        neg_csv_path = sample_dir / f"{self.config.sample}.neg.csv"

        neg_bed = pd.DataFrame(neg_bed_rows, columns=["chrom", "start", "end", "length", "id"])
        neg_csv = pd.DataFrame(neg_csv_rows, columns=["id", "fragN", "region", "length", "seq"])

        neg_bed.to_csv(neg_bed_path, sep="\t", header=False, index=False)
        neg_csv.to_csv(neg_csv_path, sep="\t", index=False)

        logger.info(f"  转换完成: {len(csv_rows)} 个 eccDNA, {len(neg_csv_rows)} 个线性 DNA")
        logger.info(f"  pos.bed: {pos_bed_path}")
        logger.info(f"  neg.bed: {neg_bed_path}")

    def _convert_region_to_lite_format(self, eccdna_fasta: Path) -> None:
        """兼容旧接口：sim-region -> lite 格式"""
        self._convert_region_format(eccdna_fasta)

    def _run_readsim(self, eccdna_fasta: Path) -> None:
        """运行读段模拟（使用 sim-region 的 eccDNA）"""
        logger.info("-" * 40)
        logger.info("步骤 2: 读段模拟 (readsim)")
        logger.info("-" * 40)

        from .readsim import libsim, fqsim

        readsim_cfg = self.config.readsim.params
        platforms = self.config.readsim.platforms

        # 将 sim-region 输出转换为读段模拟格式
        logger.info("使用 sim-region 生成的 eccDNA...")
        self._convert_region_format(eccdna_fasta)

        # 运行 libsim
        logger.info("运行 libsim...")
        libsim(
            sample=self.config.sample,
            reference=self.reference,
            path=str(self.sequencing_dir),
            seed=self.config.seed,
            meancov=readsim_cfg.meancov,
            amp=readsim_cfg.amp,
            min_repeats=readsim_cfg.min_repeats,
            threads=self.config.threads,
            rca_output_dir=str(self.rca_dir),
        )

        # 运行 fqsim
        logger.info("运行 fqsim...")
        csv_path = self.rca_dir / f"{self.config.sample}.lib.csv"

        fqsim(
            sample=self.config.sample,
            csv=str(csv_path),
            path=str(self.sequencing_dir),
            seed=self.config.seed,
            thread=self.config.threads,
            skip_sr=self.skip_sr or not platforms.ngs,
            skip_hifi=self.skip_hifi or not platforms.hifi,
            skip_ont=self.skip_ont or not platforms.ont,
            ont_mean=readsim_cfg.ont.mean_length,
            ont_std=readsim_cfg.ont.std_length,
            ont_model=readsim_cfg.ont.model,
            hifi_sample_fastq=readsim_cfg.hifi.sample_fastq,
            hifi_mode=readsim_cfg.hifi.mode,
            hifi_profile_id=readsim_cfg.hifi.profile_id,
            hifi_profile_root=readsim_cfg.hifi.profile_root,
            hifi_len_min=readsim_cfg.hifi.len_min,
            hifi_len_peak_min=readsim_cfg.hifi.len_peak_min,
            hifi_len_peak_max=readsim_cfg.hifi.len_peak_max,
            hifi_len_max=readsim_cfg.hifi.len_max,
            hifi_qmin=readsim_cfg.hifi.qmin,
            hifi_qmean=readsim_cfg.hifi.qmean,
            hifi_qsd=readsim_cfg.hifi.qsd,
            sr_mean=readsim_cfg.ngs.insert_mean,
            sr_std=readsim_cfg.ngs.insert_std,
            sr_readlen=readsim_cfg.ngs.read_length,
            sr_platform=readsim_cfg.ngs.platform,
            generate_truth=True,
        )

        logger.info("读段模拟完成")

    def _run_readsim_multi_coverage(self, eccdna_fasta: Path) -> None:
        """运行多覆盖度读段模拟"""
        from .readsim import libsim, fqsim

        readsim_cfg = self.config.readsim.params
        platforms = self.config.readsim.platforms

        # 单覆盖度时直接调用原方法
        if len(self.coverage_list) == 1:
            self._run_readsim(eccdna_fasta)
            return

        # 多覆盖度模式
        logger.info("-" * 40)
        logger.info(f"步骤 2: 多覆盖度读段模拟 ({len(self.coverage_list)} 个覆盖度)")
        logger.info(f"覆盖度列表: {', '.join([f'{self._fmt_cov(c)}X' for c in self.coverage_list])}")
        logger.info("-" * 40)

        # RCA 扩增只运行一次（使用第一个覆盖度作为基准）
        base_cov = self.coverage_list[0]
        base_cov_str = self._fmt_cov(base_cov)

        # 使用第一个覆盖度目录进行格式转换
        first_cov_dir = self.output_dir / f"sequencing_{base_cov_str}X"
        self.sequencing_dir = first_cov_dir

        logger.info("")
        logger.info("运行 RCA 扩增（只运行一次）...")
        self._convert_region_format(eccdna_fasta)

        libsim(
            sample=self.config.sample,
            reference=self.reference,
            path=str(first_cov_dir),
            seed=self.config.seed,
            meancov=base_cov,
            amp=readsim_cfg.amp,
            min_repeats=readsim_cfg.min_repeats,
            threads=self.config.threads,
            rca_output_dir=str(self.rca_dir),
        )

        # 读取基准 RCA 数据
        import pandas as pd
        base_csv_path = self.rca_dir / f"{self.config.sample}.lib.csv"
        base_lib_df = pd.read_csv(base_csv_path, sep='\t')

        # 验证必要的列存在
        required_cols = ['realcov', 'tempcov']
        missing_cols = [col for col in required_cols if col not in base_lib_df.columns]
        if missing_cols:
            raise ValueError(f"lib.csv 缺少必要的列: {missing_cols}")

        for i, cov in enumerate(self.coverage_list, 1):
            cov_str = self._fmt_cov(cov)
            cov_dir = self.output_dir / f"sequencing_{cov_str}X"
            self.sequencing_dir = cov_dir

            logger.info("")
            logger.info(f"[{i}/{len(self.coverage_list)}] 处理覆盖度 {cov_str}X")
            logger.info(f"输出目录: {cov_dir}")

            # 根据覆盖度比例调整 tempcov，生成该覆盖度的 .lib.csv
            cov_scale = cov / base_cov
            cov_lib_df = base_lib_df.copy()
            cov_lib_df['realcov'] = cov_lib_df['realcov'] * cov_scale
            cov_lib_df['tempcov'] = cov_lib_df['tempcov'] * cov_scale
            cov_csv_path = cov_dir / f"{self.config.sample}.lib.csv"
            cov_lib_df.to_csv(cov_csv_path, index=None, sep='\t')
            logger.info(f"生成覆盖度调整后的 lib.csv: {cov_csv_path}")

            # 运行 fqsim
            logger.info(f"运行 fqsim (覆盖度: {cov_str}X)...")
            csv_path = cov_csv_path

            fqsim(
                sample=self.config.sample,
                csv=str(csv_path),
                path=str(cov_dir),
                seed=self.config.seed,
                thread=self.config.threads,
                skip_sr=self.skip_sr or not platforms.ngs,
                skip_hifi=self.skip_hifi or not platforms.hifi,
                skip_ont=self.skip_ont or not platforms.ont,
                ont_mean=readsim_cfg.ont.mean_length,
                ont_std=readsim_cfg.ont.std_length,
                ont_model=readsim_cfg.ont.model,
                hifi_sample_fastq=readsim_cfg.hifi.sample_fastq,
                hifi_mode=readsim_cfg.hifi.mode,
                hifi_profile_id=readsim_cfg.hifi.profile_id,
                hifi_profile_root=readsim_cfg.hifi.profile_root,
                hifi_len_min=readsim_cfg.hifi.len_min,
                hifi_len_peak_min=readsim_cfg.hifi.len_peak_min,
                hifi_len_peak_max=readsim_cfg.hifi.len_peak_max,
                hifi_len_max=readsim_cfg.hifi.len_max,
                hifi_qmin=readsim_cfg.hifi.qmin,
                hifi_qmean=readsim_cfg.hifi.qmean,
                hifi_qsd=readsim_cfg.hifi.qsd,
                sr_mean=readsim_cfg.ngs.insert_mean,
                sr_std=readsim_cfg.ngs.insert_std,
                sr_readlen=readsim_cfg.ngs.read_length,
                sr_platform=readsim_cfg.ngs.platform,
                generate_truth=True,
            )

            logger.info(f"覆盖度 {cov_str}X 完成")

        logger.info("")
        logger.info("所有覆盖度读段模拟完成")


def show_simulation_plan(
    config: UnifiedSimulateConfig,
    output_dir: str,
    reference: str,
    skip_region: bool,
    input_eccdna: Optional[str],
    coverage_list: Optional[List[float]] = None,
) -> None:
    """显示模拟执行计划（dry-run 模式）"""
    coverage_list = coverage_list or [config.readsim.params.meancov]

    print("\n" + "=" * 60)
    print("eccDNA 模拟执行计划 (dry-run)")
    print("=" * 60)

    print(f"\n样本名称: {config.sample}")
    print(f"输出目录: {output_dir}")
    print(f"参考基因组: {reference}")
    print(f"线程数: {config.threads}")
    print(f"随机种子: {config.seed}")

    print("\n" + "-" * 40)
    print("步骤 1: 区域模拟 (sim-region)")
    print("-" * 40)
    if skip_region:
        print(f"  [跳过] 使用已有 eccDNA: {input_eccdna}")
    else:
        region = config.region
        print(f"  UeccDNA 数量: {region.num_unique}")
        print(f"  MeccDNA 数量: {region.num_multi}")
        print(f"  CeccDNA 数量: {region.num_chimeric}")
        print(f"  长度分布峰值: {region.length_distribution.mode} bp")

    print("\n" + "-" * 40)
    print("步骤 2: 读段模拟")
    print("-" * 40)
    print(f"  启用: {'是' if config.readsim.enabled else '否'}")

    platforms = config.readsim.platforms
    print(f"  平台: NGS={platforms.ngs}, HiFi={platforms.hifi}, ONT={platforms.ont}")

    if config.readsim.enabled:
        params = config.readsim.params
        print(f"\n  [读段模拟]")
        print(f"    使用 sim-region 生成的 eccDNA")
        print(f"    RCA 扩增长度: {params.amp} bp")

        # 多覆盖度支持
        def _fmt_cov(c: float) -> str:
            """格式化覆盖度：整数显示为整数，小数保留一位"""
            return f"{int(c)}" if c == int(c) else f"{c:.1f}"

        if len(coverage_list) == 1:
            print(f"    平均覆盖度: {_fmt_cov(coverage_list[0])}X")
        else:
            print(f"    多覆盖度模式: {len(coverage_list)} 个覆盖度")
            cov_str = ", ".join([f"{_fmt_cov(c)}X" for c in coverage_list])
            print(f"    覆盖度列表: {cov_str}")
            print(f"    输出目录:")
            for cov in coverage_list:
                print(f"      - {output_dir}/sequencing_{_fmt_cov(cov)}X/")
    else:
        print(f"\n  [读段模拟] 已跳过 (--skip-readsim)")

    # 验证配置
    warnings = config.validate()
    if warnings:
        print("\n" + "-" * 40)
        print("配置警告:")
        print("-" * 40)
        for w in warnings:
            print(f"  - {w}")

    print("\n" + "=" * 60)
    print("执行计划显示完成 (未实际运行)")
    print("=" * 60 + "\n")
