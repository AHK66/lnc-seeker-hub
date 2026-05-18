# lnc-seeker-collect: Algorithmic Process

This document describes the step-by-step algorithmic process performed by `lnc-seeker-collect` to extract, remap, and organize genomic reads from multiple BAM sources into region-specific collections.

## 1. Configuration and Environment Setup
The process begins by parsing the input configuration (`config.cfg`) and any command-line overrides.
- **Gene Resolution**: It identifies all gene names requested across all `gene_region` entries.
- **GTF Querying**: It queries the provided GTF file to find the genomic coordinates (chromosome, start, end) for every unique gene.
- **Region Computation**: 
  - For `gene_region` entries, it calculates the "overall region" spanned by all listed genes plus the specified offset.
  - For `region` entries, it parses the standard `chr:start-end` format.
  - All chromosome names are translated using the `assembly_report` to ensure consistency with the input BAM headers.

## 2. Extraction Tag Generation
To prevent filename conflicts and allow a single global dictionary, the tool generates unique identifiers for every requested region:
- **Descriptor Hash**: For each region entry, it creates a string descriptor comprising the gene names (or region string), the offset, and the resolved coordinates.
- **Stable Hashing**: It computes a deterministic 8-character hex hash of this descriptor.
- **Tag Construction**: A human-readable prefix (the first gene name) is combined with the hash (e.g., `FAM72A_a1b2c3d4`).

## 3. Distributed Read Collection
The tool iterates through every defined **Cohort**. For each cohort:
1. **Header Synchronization**: It reads the header of the first BAM file in the cohort to establish the "Output Reference Space."
2. **Parallel/Batch Processing**: It opens every BAM file associated with the cohort.
3. **Targeted Querying**: For each open BAM, it performs an indexed query for every requested region.
4. **Read Processing**:
   - **Filtering**: Skips secondary alignments (if configured) and reads with no alignment position.
   - **Boundary Check**: Verifies if the read actually falls within the requested coordinate window.
   - **Reference Remapping**: Translates the Reference Sequence ID (RID) from the source BAM's header space into the output BAM's header space.
   - **Endpoint Annotation**: If enabled, attaches the `EP` tag (`L`, `R`, or `LR`) to indicate which read ends fall within the collection window.
5. **Per-Region Deduplication**: 
   - Reads are placed into "buckets" specific to each region.
   - Within each bucket, reads are deduplicated based on `(Name, Start, Flags)`.
   - *Note: A read overlapping two requested regions will appear in both buckets.*

### Read Inclusion and Overlap Logic
To ensure complete genomic context, the collection algorithm implements an inclusive overlap policy:
- **Inclusive Selection**: A read is collected if any portion of its alignment overlaps the requested region. A read spanning `[s, e]` is included when `s <= region_end` and `e >= region_start`.
- **Preservation of Spanning Reads**: This approach intentionally captures reads that start outside the region but end inside (or vice versa), as well as reads that completely span the region. This is essential for visualizing long-range splice junctions and maintaining accurate coverage at the region boundaries.
- **Full Context Retention**: If a read meets the inclusion criteria, the entire record is preserved as-is. This allows downstream tools to visualize the full extent of the read even beyond the extraction coordinate boundaries.

## 4. BAM Generation and Indexing
For every Cohort-Region pair that contains reads:
- **Streaming Sort**: The collected reads are sorted by their remapped coordinate position to satisfy BAM indexing requirements.
- **BAM Writing**: A new BAM file is written to the output directory using the naming scheme: `{CohortName}_{RegionTag}.bam`.
- **BAI Indexing**: A companion `.bai` index is generated immediately for fast downstream access.

## 5. Global Manifest Generation
Finally, the tool generates a single `dictionary.json` in the output directory:
- **Hierarchical Indexing**: The JSON is keyed by the requested region string/gene names.
- **Metadata Storage**: Each section contains the specific genomic region and offset used for that collection.
- **Path Mapping**: Each section contains a `cohorts` map, where cohort names point to the path of their respective region-specific data files.

### 6. Metrics and Normalization
To facilitate cross-cohort comparison, the `dictionary.json` includes quantitative metrics for each region:

- **num_reads**: The count of unique primary/secondary alignments (after filtering) captured for the region.
- **avg_coverage_per_sample**: A normalized measure of data density across the region.

#### Coverage Formula
The average coverage per sample is calculated using the total base-length of all aligned segments (exons) within the cache:

$$Coverage = \frac{\sum_{i=1}^{R} (\sum_{j=1}^{S_i} L_{i,j})}{L_{region} \times N_{samples}}$$

Where:
- $R$ is the total number of reads.
- $S_i$ is the number of segments in read $i$.
- $L_{i,j}$ is the length of segment $j$ in read $i$.
- $L_{region}$ is the span of the requested extraction region ($End - Start + 1$).
- $N_{samples}$ is the number of BAM files (samples) defined for the cohort.

This metric allows users to identify regions with significantly different data depth across tissues or conditions, regardless of the number of samples processed.

This architecture ensures that one runs of the tool can produce granular, collision-free data sets that are immediately ready for visualization or further analysis.
