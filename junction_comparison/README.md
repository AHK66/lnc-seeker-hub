# Junction Comparison: lncRNA Seeker Hub vs StringTie

This directory contains tools to compare cohort-level splice-junction evidence
detected by lncRNA Seeker Hub against junctions present in transcript models
assembled by StringTie from the same RNA-seq BAM files.

## Workflow

1. Create a config file (see `config_comparison_example.cfg`)
2. Run StringTie on each BAM in the cohort independently (de novo, no `-G`)
3. Merge StringTie results with `stringtie --merge`
4. Extract junctions from the merged StringTie GTF
5. Resolve gene region names to genomic coordinates from the reference GTF
6. Run lnc-seeker per-region analysis on the same BAM files
7. Compare junctions at multiple read-count thresholds (1, 2, 5, 10; stranded + unstranded)
8. Classify junctions against RefSeq and GENCODE lncRNA annotations
9. Produce an enriched aggregated summary

## Prerequisites

- Python 3.10+
- StringTie 3.0+ (path specified in config)
- lnc_seeker Python module (`maturin develop --release`)
- Python packages: `pip install -r requirements.txt`

## Config File

Create a `config_comparison.cfg` in your analysis directory.
See `config_comparison_example.cfg` for a complete example.

**Important:** StringTie 3.0 cannot parse gzipped RefSeq GTFs (reads them as binary).
lnc-seeker requires BGZF-compressed GTF. Use separate entries:
- `gtf` — compressed (`.gtf.gz`) for lnc-seeker
- `gtf_stringtie` — uncompressed (`.gtf`) for StringTie / gene region resolution
If omitted, `gtf_stringtie` defaults to `gtf`.

Inline comments on `value` lines must be avoided in the config file.
Use `#` lines (before each key=value entry) for comments instead.

## Usage

### Full pipeline

```bash
./run_comparison.sh /path/to/config_comparison.cfg /path/to/output_dir
```

### Manual per-step workflow

```bash
# 1. Run StringTie on each BAM (de novo, no -G flag)
for bam in sample1.bam sample2.bam; do
    sample=$(basename "$bam" .bam)
    stringtie -o "stringtie_${sample}.gtf" "$bam"
done

# 2. Merge StringTie results
stringtie --merge -o merged.gtf stringtie_*.gtf

# 3. Extract junctions from StringTie GTF
python3 extract_stringtie_junctions.py merged.gtf -o stringtie_junctions.tsv

# 4. Resolve gene regions from GTF
python3 parse_config.py config_comparison.cfg --resolve-genes

# 5. Run lnc-seeker and export aggregated junctions
python3 export_lncseeker_junctions.py config_comparison.cfg \
    --region "chr11:65497620-65501529" \
    -o lncseeker_junctions.tsv

# 6. Compare (region-filtered)
python3 compare_junctions.py \
    --lncseeker lncseeker_junctions.tsv \
    --stringtie stringtie_junctions.tsv \
    --region "chr11:65497620-65501529" \
    -o comparison_results.csv

# 7. Classify against known annotations
python3 extract_gtf_junctions.py refseq.gtf.gz -o refseq_junctions.tsv
python3 extract_gtf_junctions.py gencode_lncrna.gtf.gz -o gencode_junctions.tsv
python3 classify_junctions.py \
    --known-refseq refseq_junctions.tsv \
    --known-gencode gencode_junctions.tsv \
    --lncseeker lncseeker_junctions.tsv \
    --stringtie stringtie_junctions.tsv \
    --region "chr11:65497620-65501529" \
    -o classification.csv

```

## Output Structure

```
output_dir/
├── stringtie/
│   ├── <sample>.gtf            # per-sample StringTie assemblies
│   ├── merged.gtf              # merged transcript assembly
│   ├── junctions.tsv           # extracted junctions from merged GTF
│   ├── refseq_junctions.tsv    # RefSeq annotation junctions
│   └── gencode_junctions.tsv   # GENCODE lncRNA annotation junctions
├── lncseeker/
│   ├── per_region/
│   │   └── <gene>_junctions.tsv   # per-gene aggregated junctions
│   └── <gene>_junctions.tsv
├── comparison/
│   ├── <gene>/                     # per-region results
│   │   └── results.csv                    # overlap stats at each threshold
│   ├── classification/             # per-region annotation classification
│   │   ├── <gene>.csv
│   │   └── ...
│   └── aggregated_summary.csv      # enriched summary with annotation fields
├── config_resolved.json            # resolved config with gene coordinates
└── logs/
```

## Coordinate System

All junctions are in 0-based half-open format: `(chrom, start, end, strand)`.
StringTie's 1-based coordinates are converted during extraction: `end_0based = jend - 1`.
Chromosome names are normalized across sources (RefSeq `NC_` accession formats
are mapped to the corresponding chromosome numbers).
