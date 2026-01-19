"""PacBio HiFi read generation for RCA simulation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple
import logging

import numpy as np

from ..config import SimConfig
from ..debranch import LinearMoleculePool
from ..models import SequencedRead, SegmentType, compute_junction_covered_possible

logger = logging.getLogger(__name__)


@dataclass
class HiFiLibraryStats:
    """Summary stats for HiFi read generation."""

    total_reads: int = 0
    mean_length: float = 0.0
    min_length: int = 0
    max_length: int = 0
    background_reads: int = 0
    chimera_reads: int = 0

    def summary(self) -> str:
        if self.total_reads <= 0:
            return "HiFi reads: 0"
        return (
            f"HiFi reads: {self.total_reads}, "
            f"len={self.min_length}-{self.max_length} (mean={self.mean_length:.1f}), "
            f"background={self.background_reads}, chimera={self.chimera_reads}"
        )


class HiFiLibraryGenerator:
    """Generate HiFi reads directly from linear molecules."""

    def __init__(self, config: SimConfig, rng: Optional[np.random.Generator] = None):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._read_counter = 0

    def _sample_read_length(self, mol_len: int) -> int:
        const = self.config.constants
        min_len = max(1, int(const.hifi_min_length))
        max_len = max(min_len, int(const.hifi_max_length))
        mean_len = float(const.hifi_target_size_mean)
        std_len = float(const.hifi_target_size_std)

        if mol_len <= 0:
            return 0

        effective_max = min(max_len, mol_len)
        effective_min = min(min_len, effective_max)

        if std_len <= 0:
            candidate = int(round(mean_len))
            return max(effective_min, min(effective_max, candidate))

        for _ in range(100):
            candidate = int(round(self.rng.normal(mean_len, std_len)))
            if effective_min <= candidate <= effective_max:
                return candidate

        candidate = int(round(self.rng.normal(mean_len, std_len)))
        return max(effective_min, min(effective_max, candidate))

    def _build_read(
        self,
        pool: LinearMoleculePool,
        mol,
        read_len: int
    ) -> Optional[SequencedRead]:
        mol_len = mol.total_length
        if mol_len <= 0 or read_len <= 0:
            return None

        max_start = mol_len - read_len
        start = 0 if max_start <= 0 else int(self.rng.integers(0, max_start + 1))

        seq = pool.get_subsequence(mol, start, read_len)
        segments = mol.extract_subregion(start, read_len)
        if not segments:
            return None

        is_background = mol.is_background or all(
            seg.segment_type == SegmentType.BACKGROUND for seg in segments
        )

        source_ecc_ids: List[str] = []
        repeat_count_by_source: dict = {}
        if not is_background:
            for seg in segments:
                if seg.segment_type == SegmentType.BACKGROUND:
                    continue
                if seg.ecc_id not in source_ecc_ids:
                    source_ecc_ids.append(seg.ecc_id)
                ecc = pool.ecc_db.get(seg.ecc_id)
                ecc_len = ecc.length if ecc is not None else 0
                if ecc_len > 0:
                    repeat_count_by_source[seg.ecc_id] = (
                        repeat_count_by_source.get(seg.ecc_id, 0.0) + seg.length / ecc_len
                    )

        repeat_count_truth = sum(repeat_count_by_source.values()) if repeat_count_by_source else 0.0

        chimera_breakpoints: List[int] = []
        pos = 0
        for seg in segments:
            if seg.segment_type == SegmentType.CHIMERA:
                chimera_breakpoints.append(pos)
            pos += seg.length

        has_inter_chimera = bool(chimera_breakpoints) or len(source_ecc_ids) > 1

        read_index = self._read_counter
        self._read_counter += 1
        read_id = f"{mol.molecule_id}_hifi_{read_index}"

        return SequencedRead(
            read_id=read_id,
            platform="HiFi",
            sequence=seq,
            quality="I" * len(seq),
            source_molecule_id=mol.molecule_id,
            source_ecc_ids=source_ecc_ids,
            segments=segments,
            is_background=is_background,
            background_chrom=mol.background_chrom if is_background else None,
            background_start=mol.background_start if is_background else None,
            background_end=mol.background_end if is_background else None,
            repeat_count_truth=repeat_count_truth,
            repeat_count_by_source=repeat_count_by_source,
            has_inter_chimera=has_inter_chimera,
            chimera_breakpoints=chimera_breakpoints,
            junction_covered_possible=compute_junction_covered_possible(segments),
        )

    def generate(
        self,
        pool: LinearMoleculePool,
        target_reads: int
    ) -> Tuple[List[SequencedRead], HiFiLibraryStats]:
        if target_reads <= 0:
            return [], HiFiLibraryStats()
        if pool.total_molecules <= 0:
            raise ValueError("Cannot generate HiFi reads from empty molecule pool")

        molecules = pool.sample(target_reads, self.rng)
        reads: List[SequencedRead] = []
        lengths: List[int] = []

        for mol in molecules:
            read_len = self._sample_read_length(mol.total_length)
            read = self._build_read(pool, mol, read_len)
            if read is None:
                continue
            reads.append(read)
            lengths.append(len(read.sequence))

        stats = HiFiLibraryStats(total_reads=len(reads))
        if lengths:
            stats.min_length = int(min(lengths))
            stats.max_length = int(max(lengths))
            stats.mean_length = float(np.mean(lengths))
            stats.background_reads = sum(1 for r in reads if r.is_background)
            stats.chimera_reads = sum(1 for r in reads if r.has_inter_chimera)

        return reads, stats
