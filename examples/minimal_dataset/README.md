# Minimal Example Dataset

A small, fully public dataset for running **lncRNA Seeker Hub** end to end and
confirming that your installation is configured and executed correctly. It reuses
a three-sample subset of the public **GSE235167** patient-derived-xenograft
breast-cancer RNA-seq cohort, restricted to the six loci analysed in
**Supplementary Note 6**.

The main analysis relies on restricted-access TCGA data. This example is built
from a public GEO accession so that any user can pull it straight from the
repository and run the tool without any data licensing.

## Provenance

- **Source** — GSE235167 (PDX breast cancer, Illumina HiSeq 4000, ~300 bp paired-end).
- **Samples** — three of the 20 PDX samples used in Supplementary Note 6, chosen
  to span the two molecular subtypes and primary vs metastatic tissue sites:
  `SRR24952131` (Basal-like TNBC, mammary tumor, untreated),
  `SRR24952181` (ER+, mammary tumor, estrogen),
  `SRR24952121` (Basal-like TNBC, lung metastasis, untreated).
- **Loci (±500 bp padding)** — `ACTB`, `PCA3`, `CD44`, `MALAT1`, `GAPDH`, `TP53`
  (the same six loci from Supplementary Note 6).
- **BAMs** contain *only* the reads mapping to these six loci (sub-sampled from
  the full genome), which is why they are a few tens of MB instead of gigabytes.
- **Annotation** — a RefSeq GTF subset containing exactly these six genes
  (BGZF-compressed + tabix-indexed, as the tool requires).
- **`assembly_report.txt`** — the RefSeq GRCh38.p14 assembly report, used to map
  `NC_` accession names (GTF) to `chr` names (BAM).

## Contents

| File | Description | Size |
|---|---|---|
| `bam/SRR24952131_…out.bam` (+`.bai`) | Basal TNBC, mammary, 6 loci | ~22 MB |
| `bam/SRR24952181_…out.bam` (+`.bai`) | ER+, mammary, 6 loci | ~27 MB |
| `bam/SRR24952121_…out.bam` (+`.bai`) | Basal TNBC, lung met, 6 loci | ~20 MB |
| `annotation.gtf.gz` (+`.tbi`) | RefSeq GTF subset: the six genes | ~32 KB |
| `assembly_report.txt` | GRCh38.p14 assembly report | ~79 KB |
| `config_example.cfg` | Ready-to-run `lnc-seeker-collect` config | — |
| `config.json.example` | Optional: app config to view results | — |
| `run_example.sh` | End-to-end smoke test | — |

Total: **~68 MB** (largest single file ~27 MB, well under GitHub's 100 MB limit).

## Requirements

- The `lnc-seeker-collect` binary. Build it from the repository root:
  `cargo build -p lnc-seeker-collect`.
- **No reference genome is needed.** Mismatch indexing is optional, so the
  ~1 GB GRCh38 FASTA is not required to run this example.

## Quick start (smoke test)

```bash
cd examples/minimal_dataset
./run_example.sh
```

This runs `lnc-seeker-collect` against `config_example.cfg`, builds one compact
cache per locus in `collect_output/`, and verifies the outputs. A successful run
prints `SUCCESS: minimal dataset ran end to end.` and produces:

- `collect_output/dictionary.json` — per-locus metadata (reads, coverage, paths).
- `collect_output/GSE235167_<gene>_<hash>.lnc_cache.bin` — one cache per locus
  (six files).

**Built-in sanity check.** `PCA3` is a prostate-specific lncRNA, so in these
breast samples its cache is tiny (a few dozen reads) compared with `GAPDH`. The
script prints both sizes; if you ever see a large `PCA3` cache in breast data,
the inputs are misconfigured.

## Open the results in the interactive app (optional)

After `run_example.sh`, you can visualise the collected data with the Bokeh
front end. From this directory:

```bash
cp config.json.example config.json          # config.json is git-ignored
python -m bokeh serve --show ../../lnc_seeker_server.py --session-token-expiration 3600
```

`config.json` is read from the working directory; the app reads the collection
from `collect_output/`. See the repository [README](../../README.md) for the
full application quickstart.

## Regenerating the dataset

For reference, the BAMs and GTF were produced from the full public files as:

```bash
# 1) Sub-sample each full BAM to the six loci (2 kb margin), then index
regions="chr7:5524648-5533101 chr9:76761936-76790069 chr11:35136671-35234902 \
         chr11:65495238-65509016 chr12:6532017-6540871 chr17:7665921-7689990"
for s in SRR24952131 SRR24952181 SRR24952121; do
    samtools view "$FULL_BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam" $regions \
        -o "bam/${s}_Aligned.sortedByCoord.out.bam"
    samtools index "bam/${s}_Aligned.sortedByCoord.out.bam"
done

# 2) Extract the six genes from the full RefSeq GTF, then BGZF + tabix
zcat full_refseq.gtf.gz | grep -E ' gene "(ACTB|PCA3|CD44|MALAT1|GAPDH|TP53)";' \
    > annotation.gtf
bgzip annotation.gtf && tabix -s 1 -b 4 -e 5 annotation.gtf.gz
```
