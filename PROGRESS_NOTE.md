# Progress Note

Date: 2024-12-30

Work done:
- Removed empty directories accidentally created by CLI flags:
  --cov-hifi, --seed, --skip-ngs, --skip-ont, -i, -o, -t, -v, 100, 4, 42.

Findings (not yet changed in code):
- `Data/truth.HiFi.tsv` does not match `Data/reads_hifi.fastq`.
  - HiFi read IDs in FASTQ are `hifi_read_*`.
  - `truth.HiFi.tsv` uses `w*_hifi_read_*` and has 100000 rows.
  - `Data/truth.tsv` matches `reads_hifi.fastq`.
- Current pipeline only writes one truth file:
  `src/ecctoolkit/simulate/reads.py` writes `truth.tsv` (or jsonl), not per-platform truth.
- `--cov-ngs 1000` produces 1000 total NGS reads (R1/R2 are 500 each),
  because NGS target pairs = target_reads // 2 in
  `src/ecctoolkit/simulate/rca_readsim/library/ngs.py`.

Open questions:
- Should we add `truth.NGS/HiFi/ONT.tsv` outputs?
- Should `--cov-ngs` represent pairs rather than total reads?
