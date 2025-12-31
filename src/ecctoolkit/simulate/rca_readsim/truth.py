"""
Truth追踪与输出模块

实现§11：truth规范
- TSV格式输出
- JSONL格式输出
- 完整的segments追踪
"""

from typing import List, Dict, Optional, TextIO, Union
from pathlib import Path
import json
import csv

from .models import SequencedRead, Segment


class TruthWriter:
    """Truth输出写入器"""
    
    def __init__(self, output_path: Union[str, Path], format: str = "tsv"):
        """
        Args:
            output_path: 输出文件路径
            format: 输出格式 ("tsv" 或 "jsonl")
        """
        self.output_path = Path(output_path)
        self.format = format.lower()
        self._file: Optional[TextIO] = None
        self._csv_writer = None
        self._header_written = False
        
        # TSV列定义
        self._columns = [
            "read_id",
            "platform",
            "source_type",  # eccDNA or background
            "source_molecule_id",
            "source_ecc_id",
            "source_ecc_ids",
            "repeat_count_truth",
            "repeat_count_by_source_json",
            "segments_json",
            "has_inter_chimera",
            "has_branch_chimera",
            "chimera_breakpoints",
            "covered_branch_anchor",
            "truncated",
            "truncation_reason",
            "junction_covered_possible",
            "is_paired",
            "mate_id",
            "read_number",
            "insert_size",
            "read_length",
            # Background DNA specific columns
            "background_chrom",
            "background_start",
            "background_end"
        ]
    
    def open(self):
        """打开输出文件"""
        self._file = open(self.output_path, 'w', newline='', encoding='utf-8')
        if self.format == "tsv":
            self._csv_writer = csv.DictWriter(
                self._file, 
                fieldnames=self._columns,
                delimiter='\t',
                extrasaction='ignore'
            )
    
    def close(self):
        """关闭输出文件"""
        if self._file:
            self._file.close()
            self._file = None
    
    def __enter__(self):
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def write_read(self, read: SequencedRead):
        """
        写入单条read的truth记录
        
        Args:
            read: SequencedRead对象
        """
        if self._file is None:
            raise RuntimeError("Writer not opened")
        
        record = self._read_to_record(read)
        
        if self.format == "jsonl":
            self._file.write(json.dumps(record, ensure_ascii=False) + "\n")
        else:  # tsv
            if not self._header_written:
                self._csv_writer.writeheader()
                self._header_written = True
            self._csv_writer.writerow(record)
    
    def write_reads(self, reads: List[SequencedRead]):
        """
        批量写入reads
        
        Args:
            reads: SequencedRead列表
        """
        for read in reads:
            self.write_read(read)
    
    def _read_to_record(self, read: SequencedRead) -> dict:
        """将read转换为记录字典"""
        # 基本信息
        record = {
            "read_id": read.read_id,
            "platform": read.platform,
            "source_type": "background" if read.is_background else "eccDNA",
            "source_molecule_id": read.source_molecule_id,
            "source_ecc_id": read.source_ecc_ids[0] if len(read.source_ecc_ids) == 1 else "multi",
            "source_ecc_ids": ",".join(read.source_ecc_ids),
            "repeat_count_truth": read.repeat_count_truth,
            "has_inter_chimera": int(read.has_inter_chimera),
            "has_branch_chimera": int(read.has_branch_chimera),
            "chimera_breakpoints": ",".join(map(str, read.chimera_breakpoints)),
            "covered_branch_anchor": int(read.covered_branch_anchor),
            "truncated": int(read.truncated),
            "truncation_reason": read.truncation_reason,
            "junction_covered_possible": int(read.junction_covered_possible),
            "is_paired": int(read.is_paired),
            "mate_id": read.mate_id or "",
            "read_number": read.read_number,
            "insert_size": read.insert_size,
            "read_length": len(read.sequence),
            # Background DNA specific fields
            "background_chrom": read.background_chrom or "",
            "background_start": read.background_start or 0,
            "background_end": read.background_end or 0
        }

        repeat_by_source = read.repeat_count_by_source or {}
        if self.format == "jsonl":
            record["repeat_count_by_source"] = repeat_by_source
        else:
            record["repeat_count_by_source_json"] = json.dumps(repeat_by_source, ensure_ascii=False)

        # Segments序列化
        segments_data = [s.to_dict() for s in read.segments]
        record["segments_json"] = json.dumps(segments_data, ensure_ascii=False)

        return record


class SegmentTsvWriter:
    """Streaming segment TSV writer per platform."""

    def __init__(self, output_dir: Union[str, Path], prefix: str = ""):
        self.output_path = Path(output_dir)
        self.prefix = prefix
        self._writers: Dict[str, TextIO] = {}

    def _get_writer(self, platform: str) -> TextIO:
        key = platform.lower()
        if key not in self._writers:
            tsv_path = self.output_path / f"{self.prefix}reads_{key}.segments.tsv"
            self._writers[key] = open(tsv_path, "w", encoding="utf-8")
        return self._writers[key]

    def write_read(self, read: SequencedRead) -> None:
        if not read.segments:
            return
        writer = self._get_writer(read.platform)
        source_type = "background" if read.is_background else "eccDNA"
        ecc_ids = ",".join(read.source_ecc_ids)
        repeat_json = json.dumps(read.repeat_count_by_source or {}, ensure_ascii=False)

        for seg in read.segments:
            if seg.length <= 0:
                continue
            if read.is_background and read.background_chrom:
                chrom = read.background_chrom
                base_start = read.background_start or 0
                start = base_start + int(seg.ecc_offset)
                end = start + int(seg.length)
            else:
                chrom = seg.ecc_id
                start = int(seg.ecc_offset)
                end = start + int(seg.length)

            name = (
                f"{read.read_id}|seg={seg.segment_type.value}"
                f"|eccs={len(read.source_ecc_ids)}|repeat={read.repeat_count_truth}"
            )
            score = int(read.repeat_count_truth)
            strand = seg.strand.value

            writer.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t"
                f"{read.platform}\t{source_type}\t{ecc_ids}\t{repeat_json}\n"
            )

    def write_reads(self, reads: List[SequencedRead]) -> None:
        for read in reads:
            self.write_read(read)

    def close(self) -> None:
        for handle in self._writers.values():
            handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class TruthReader:
    """Truth文件读取器（用于验证/分析）"""
    
    def __init__(self, input_path: Union[str, Path]):
        self.input_path = Path(input_path)
        self.format = "jsonl" if self.input_path.suffix == ".jsonl" else "tsv"
    
    def read_all(self) -> List[dict]:
        """读取所有记录"""
        records = []
        
        with open(self.input_path, 'r', encoding='utf-8') as f:
            if self.format == "jsonl":
                for line in f:
                    line = line.strip()
                    if line:
                        records.append(json.loads(line))
            else:  # tsv
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    # 解析segments
                    if 'segments_json' in row and row['segments_json']:
                        row['segments'] = json.loads(row['segments_json'])
                    if 'repeat_count_by_source_json' in row and row['repeat_count_by_source_json']:
                        row['repeat_count_by_source'] = json.loads(row['repeat_count_by_source_json'])
                    # 转换数值字段
                    for key in ['repeat_count_truth', 'has_inter_chimera',
                                'has_branch_chimera', 'covered_branch_anchor', 'truncated',
                                'junction_covered_possible', 'is_paired', 'read_number',
                                'insert_size', 'read_length']:
                        if key in row and row[key]:
                            row[key] = int(row[key])
                    records.append(row)

        for row in records:
            if 'chimera_breakpoints' in row:
                row['chimera_breakpoints'] = _parse_chimera_breakpoints(row['chimera_breakpoints'])
        
        return records
    
    def iter_records(self):
        """迭代读取记录"""
        with open(self.input_path, 'r', encoding='utf-8') as f:
            if self.format == "jsonl":
                for line in f:
                    line = line.strip()
                    if line:
                        row = json.loads(line)
                        if 'chimera_breakpoints' in row:
                            row['chimera_breakpoints'] = _parse_chimera_breakpoints(row['chimera_breakpoints'])
                        yield row
            else:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if 'segments_json' in row and row['segments_json']:
                        row['segments'] = json.loads(row['segments_json'])
                    if 'repeat_count_by_source_json' in row and row['repeat_count_by_source_json']:
                        row['repeat_count_by_source'] = json.loads(row['repeat_count_by_source_json'])
                    if 'chimera_breakpoints' in row:
                        row['chimera_breakpoints'] = _parse_chimera_breakpoints(row['chimera_breakpoints'])
                    yield row


def _parse_chimera_breakpoints(value: Union[str, List[int], None]) -> List[int]:
    """解析chimera_breakpoints字段为int列表"""
    if value is None:
        return []
    if isinstance(value, list):
        return [int(v) for v in value]
    if isinstance(value, str):
        if not value:
            return []
        return [int(v) for v in value.split(',') if v]
    return []


def write_segments_by_platform(
    reads: List[SequencedRead],
    output_dir: Union[str, Path],
    prefix: str = "",
) -> None:
    """
    Write segment TSV files per platform from read segments.

    Output files:
      - reads_ngs.segments.tsv
      - reads_hifi.segments.tsv
      - reads_ont.segments.tsv

    Columns (TSV):
      chrom, start, end, name, score, strand, platform, source_type, source_ecc_ids,
      repeat_count_by_source_json
    """
    with SegmentTsvWriter(output_dir, prefix=prefix) as writer:
        writer.write_reads(reads)


BedWriter = SegmentTsvWriter
write_bed_by_platform = write_segments_by_platform


def compute_truth_statistics(records: List[dict]) -> dict:
    """
    计算truth统计信息
    
    Args:
        records: truth记录列表
    
    Returns:
        统计字典
    """
    stats = {
        "total_reads": len(records),
        "by_platform": {},
        "chimera_rate": 0.0,
        "branch_anchor_rate": 0.0,
        "truncation_rate": 0.0,
        "repeat_count_distribution": {},
        "ecc_coverage": {}
    }
    
    if not records:
        return stats
    
    # 按平台统计
    for rec in records:
        platform = rec.get("platform", "unknown")
        if platform not in stats["by_platform"]:
            stats["by_platform"][platform] = 0
        stats["by_platform"][platform] += 1
        
        # Chimera
        if rec.get("has_inter_chimera"):
            stats["chimera_rate"] += 1
        
        # Branch anchor
        if rec.get("covered_branch_anchor"):
            stats["branch_anchor_rate"] += 1
        
        # Truncation
        if rec.get("truncated"):
            stats["truncation_rate"] += 1
        
        # Repeat count
        rc = rec.get("repeat_count_truth", 0)
        if rc not in stats["repeat_count_distribution"]:
            stats["repeat_count_distribution"][rc] = 0
        stats["repeat_count_distribution"][rc] += 1
        
        # EccDNA coverage
        ecc_ids = rec.get("source_ecc_ids", "").split(",")
        for eid in ecc_ids:
            if eid:
                if eid not in stats["ecc_coverage"]:
                    stats["ecc_coverage"][eid] = 0
                stats["ecc_coverage"][eid] += 1
    
    # 计算比例
    n = len(records)
    stats["chimera_rate"] /= n
    stats["branch_anchor_rate"] /= n
    stats["truncation_rate"] /= n
    
    return stats
