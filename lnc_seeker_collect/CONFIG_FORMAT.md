# lnc-seeker-collect Configuration Format

The configuration file for `lnc-seeker-collect` uses a simple `key=value` format with support for comments, cohort-specific sections, and modular includes.

## General Structure

- **Comments**: Any line starting with `#` is ignored.
- **Global Settings**: Defined at the top level or within included files.
- **Sections**: Defined using `[cohort:Name]`. Keys appearing after a section header apply to that specific cohort until the next section header or the end of the file.
- **Includes**: The `include=path/to/file` directive allows splitting configurations across multiple files.

---

## Global Directives

| Key | Description |
|:---|:---|
| `gtf` | **(Required)** Path to the reference GTF file (supports `.gz`). |
| `assembly_report` | **(Required)** Path to the NCBI/RefSeq assembly report mapping sequence names to accessions. |
| `genome` | Path to the reference genome FASTA file. If provided, compact storage will include mismatch information. |
| `output_dir` | Path to a directory where collected BAM files and the `dictionary.json` will be saved. |
| `include` | Path to another configuration file to include. Paths are resolved relative to the file containing the `include` directive. |
| `annotate_endpoint_tag` | (Boolean: `true`/`false`) If enabled, adds an `EP:Z:` tag to records indicating if their start/end fall within extraction regions. |
| `skip_secondary` | (Boolean: `true`/`false`) If enabled, secondary alignments (flag `0x100`) are ignored. |
| `compress_headers` | (Boolean: `true`/`false`) Enables experimental differential header compression for read names in the output cache. Significantly reduces disk footprint for localized regions like introns. |
| `header_compression_algorithm` | Algorithm for header compression when `compress_headers` is enabled: `huffman` (Delta+Huffman, default), `zstd` (Delta+Zstd), or `delta` (pure Delta without additional compression). |

---

## Extraction Regions

Regions can be specified globally and will be applied to all cohorts.

| Key | Format | Description |
|:---|:---|:---|
| `region` | `chrom:start-end` | Manually specify a genomic range to extract. |
| `gene_region` | `GENES[:offset]` | Automatically compute a region spanning one or more genes (fetched from GTF) plus an optional padding `offset` in base pairs (defaults to `0`). `GENES` can be a single gene name or a comma-separated list. |

---

## Cohort Sections

Sections are defined as `[cohort:CohortName]`. Within a section, the following keys are recognized:

| Key | Description |
|:---|:---|
| `bam` | **(Repeatable)** Path to an input BAM file. Every `bam` entry is treated as a sample belonging to this cohort. |
| `output_bam` | Override the default output path (`output_dir/CohortName.bam`) for this specific cohort. |
| `tissue` | (Metadata) Descriptive string for the tissue type (e.g., `breast`, `liver`). |
| `status` | (Metadata) Descriptive string for the sample status (e.g., `cancer`, `normal`). |

**Note**: The number of samples (`num_samples`) is automatically calculated based on the count of `bam=` entries in the section.

---

## Example Configuration

```config
# Global Paths
gtf = /data/genome/GRCh38.gtf.gz
assembly_report = /data/genome/GRCh38_assembly_report.txt
output_dir = /data/output/extracted_results

# Modular settings
include = shared_extraction_logic.cfg

# Regions of interest
gene_region = FAM72A,SRGAP2:5000

# Cohort A
[cohort:GSE235167]
tissue = breast
status = cancer
bam = /data/raw/sample1_Aligned.bam
bam = /data/raw/sample2_Aligned.bam

# Cohort B
[cohort:Healthy_Control]
tissue = breast
status = normal
bam = /data/raw/control1_Aligned.bam
```
