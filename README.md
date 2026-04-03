# plasmidsaurus_scmpra_troubleshooting

Long-read sequencing analysis for troubleshooting scmPRA PCR amplification. Aligns Plasmidsaurus reads to reference amplicons, characterizes contaminants, and quantifies R1-primer-flanked reads to diagnose PCR failures.

## Environment
Conda env used: `plasmidsaurus-align`

Environment exports are committed under `env/`:
- `env/environment.plasmidsaurus-align.yml` (portable conda spec)
- `env/environment.plasmidsaurus-align.from-history.yml` (minimal requested specs)
- `env/environment.plasmidsaurus-align.lock.txt` (`conda --explicit` lock-style spec)
- `env/environment.plasmidsaurus-align.pip.txt` (`pip freeze`)

Refresh these files after package changes:
```bash
scripts/export_env.sh plasmidsaurus-align
```

## Implemented pipeline
Driver script:
```bash
conda activate plasmidsaurus-align
scripts/run_alignment_pipeline.sh \
  --sample <sample_id> \
  --reads data/reads/<sample>.fastq.gz \
  --reference-gb data/references/<ref>.gb \
  --threads 4
```

Non-activated shell alternative:
```bash
conda run -n plasmidsaurus-align bash scripts/run_alignment_pipeline.sh \
  --sample <sample_id> \
  --reads data/reads/<sample>.fastq.gz \
  --reference-gb data/references/<ref>.gb \
  --threads 4
```

What it runs:
1. `scripts/gb_to_fasta_and_bed.py` to convert GenBank to FASTA + BED feature intervals.
2. `minimap2` index + alignment to BAM.
3. `samtools depth -aa` for base-wise coverage.
4. `bedtools coverage` for raw per-feature overlap metrics.
5. `scripts/summarize_feature_coverage.py` for feature-level summary table + mean-depth plot.
6. `scripts/plot_depth_with_features.py` for depth line plot with annotation tracks.

Outputs:
- `results/<sample>/<sample>.bam` and `.bam.bai`
- `results/<sample>/<sample>.flagstat.txt`
- `results/<sample>/<sample>.depth.tsv`
- `results/<sample>/<sample>.feature_coverage.raw.tsv`
- `results/<sample>/<sample>.feature_coverage.summary.tsv`
- `results/<sample>/<sample>.depth_with_features.png`
- `results/<sample>/<sample>.feature_mean_depth.png`

## Human RNA-cDNA explainability (transcriptome-only)
This workflow aligns each read to a human transcriptome reference and calls whether
the read is "explained as human RNA cDNA".

Definition (default):
- `best_query_cov >= 0.50`
- `best_identity >= 0.80`

This is transcriptome-only by design (not a genome alignment workflow).

### 1) Fetch/cache transcript reference
```bash
conda run -n plasmidsaurus-align bash scripts/fetch_human_transcriptome_refs.sh --release v48
```

Outputs:
- `data/human_refs/v48/human_transcripts.fa.gz`
- `data/human_refs/v48/human_transcripts.fa`
- `data/human_refs/v48/human_transcripts.mmi`
- `data/human_refs/v48/REFERENCE_INFO.tsv`

### 2) Quantify per-read explainability
```bash
conda run -n plasmidsaurus-align python scripts/quantify_human_cdna_explained.py \
  --reads data/reads/<sample>.fastq \
  --ref-mmi data/human_refs/v48/human_transcripts.mmi \
  --threads 8 \
  --min-query-cov 0.50 \
  --min-identity 0.80 \
  --out-prefix results/<sample>/human_cdna_explainability
```

Primary outputs:
- `results/.../human_cdna_explainability.summary.tsv`
- `results/.../human_cdna_explainability.length_binned.tsv`
- `results/.../human_cdna_explainability.top_targets.tsv`
- `results/.../human_cdna_explainability.threshold_grid.tsv`
- `results/.../human_cdna_explainability.identity_hist.tsv`
- `results/.../human_cdna_explainability.query_cov_hist.tsv`
- `results/.../human_cdna_explainability.alignment.paf.gz`
- optional detailed table: `results/.../human_cdna_explainability.per_read.tsv`

Required columns in `per_read.tsv`:
- `read_id`
- `read_length`
- `best_target`
- `best_query_cov`
- `best_identity`
- `explained`

Key metrics in `summary.tsv`:
- `reads_explained_fraction = explained_reads / total_reads`
- `bases_explained_fraction = explained_bases / total_bases`
- `mapped_to_human_fraction = mapped_to_human_reads / total_reads`

## Suggested layout
```text
data/
  reads/
  references/
scripts/
results/
```

## Sample mapping table
Use this table to track which read file is mapped to which reference and whether the target is an amplicon or plasmid.

| sample_id | reads_file | reference_genbank | reference_seq_name | status (amplicon/plasmid) | notes |
|---|---|---|---|---|---|
| LFFHG8_1_pcr1_sub5 | `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq` | `data/references/gblock_f1r1.gb` | `gblock_f1r1` | amplicon | 4/4586 reads mapped (0.09%); dominated by R1-flanked rDNA contaminant (~21.5%) |
| LFFHG8_1_pcr1_sub5_mock | `data/reads/LFFHG8_1_pcr1_sub5.mock.fastq` | `data/references/gblock_f1r1.gb` | `gblock_f1r1` | amplicon | Historical smoke-test input generated from consensus FASTA |
| 3LV7H5_1_sub5 | `data/reads/3LV7H5_1_5.fastq` | -- | -- | amplicon | R1-flank + contaminant analysis only; 49.5% R1-flanked reads |
| 3LV7H5_2_sub6 | `data/reads/3LV7H5_2_6.fastq` | -- | -- | amplicon | R1-flank + contaminant analysis only; 56.6% R1-flanked reads |
