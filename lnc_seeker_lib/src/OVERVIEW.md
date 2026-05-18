# lnc_seeker_lib Source Map

> **NOTICE FOR AI AGENTS:** This file provides a high-level orientation of the library's architecture. If you perform significant updates, refactorings, or add new modules, you **must** update this file to reflect the changes.

## Core Architecture

The library is structured as a modular bioinformatics tool that processes BAM and GTF files to identify novel transcripts (lncRNAs) and visualize their coverage. It provides a C-API (via PyO3) for integration with the Python-based Bokeh frontend.

### Entry Points

- **[lib.rs](lib.rs)**: The library entry point. It declares the module structure, re-exports essential symbols, and contains the `#[pymodule]` definition and thin `#[pyfunction]` wrappers that bridge Rust logic to Python.
- **[main.rs](main.rs)**: The Command Line Interface (CLI) entry point. It uses `clap` to parse arguments and invokes `analysis::run_analysis`.

### Logic Modules

- **[analysis.rs](analysis.rs)**: Orchestrates the high-level workflow. It manages the parallel processing of BAM files via `rayon`, handles directory cleanup, and coordinates the interaction between coverage calculation and annotation lookup. It supports **BAM-free operation** by prioritizing compact-storage caches.
- **[annotations.rs](annotations.rs)**: Contains the logic for querying GTF files. It supports both high-performance indexed lookups (CSI/Tabix) via `noodles` and a streaming fallback parser.
- **[assembly.rs](assembly.rs)**: Handles `AssemblyReport` parsing. This is used to map between different chromosome naming conventions (e.g., INSDC vs. RefSeq vs. UCSC).
- **[compression.rs](compression.rs)**: High-performance Differential Huffman compression engine for read names. Implements prefix-biased delta encoding and $O(1)$ table-driven canonical decoding.
- **[coverage.rs](coverage.rs)**: The core genomics engine. It tracks depth of coverage across references and identifies junction points (splice sites) by parsing CIGAR strings. It is **metadata-agnostic**, working with cached reference lengths instead of requiring raw BAM headers.
- **[genome.rs](genome.rs)**: Provides indexed FASTA access for comparing read sequences against the reference genome to identify mismatches.
- **[processing.rs](processing.rs)**: Logic-heavy functions used primarily by the UI, including downsampling large coverage datasets for plotting and extracting specific reads that cross junctions of interest. It leverages the global `ReadProvider` for rapid access to cached data.
- **[reads_manager.rs](reads_manager.rs)**: The brain of the caching system. Manages the global `ReadProvider`, handles `LNC1` multi-block file persistence, and implements parallel disk-to-RAM hydration logic.

### Data & Configuration

- **[config.rs](config.rs)**: Defines the `Config` and `DataSelection` structures used to deserialize `config.json`. Includes custom logic for handling flexible path inputs.
- **[models.rs](models.rs)**: Centralized definition of shared data structures used across the library, such as `SampleResult`, `Annotation`, `JunctionSpan`, and `ReadInfo`.
- **[progress.rs](progress.rs)**: Implements the `SessionProgress` type and high-level progress tracking logic using atomic variables to provide real-time updates to the Python UI.

### Utilities

- **[utils.rs](utils.rs)**: Common helper functions, including GTF attribute parsing (`parse_attributes`) and interval intersection logic (`intersects`).

---

## Key Features & Logic

### Genome-Based Mismatch Discovery
Standard SAM/BAM "Match" CIGAR operations (`M`) do not distinguish between sequence identity and mismatches. `lnc_seeker_lib` implements a reliable mismatch discovery engine:
1. If a `genome_path` (indexed FASTA) is configured, the `GenomeProvider` loads the relevant reference sequences.
2. During BAM parsing, the read sequence is compared base-by-base against the reference.
3. Mismatched bases are recorded in `ReadInfo` and persisted in the `BamCache`, allowing accurate visualization of SNPs and variants in the UI.

#### Visual Representation (Full Read Layouts)
The UI provides an efficient visualization of sequence variations:
- **Consistent Styling**: Aligned segments use base colors (Blue for single, Dark Red for mates) and uniform thickness.
- **Variation Highlighting**: Segments containing mismatches or insertions are automatically displayed in a **darker shade** (e.g., Midnight Blue or Dark Red) to provide subtle but effective visual grouping.
- **High-Contrast Labels**: Sequence variations are highlighted using black boxes with white text for maximum visibility.
- **Labels**: Segments display their exact variation counts (e.g., `2X 1I`) where **X** represents mismatches and **I** represents insertions.

### Junction Quality Filtering
Splice junctions are validated by checking the "cleanliness" of their flanking segments. 
- In **[coverage.rs](coverage.rs)**, the `JunctionStore` maintains separate counts for "total" and "clean" reads. A read is clean if its flanking segments have zero mismatches and insertions (as determined by the `GenomeProvider`).
- In **[processing.rs](processing.rs)**, the `get_junction_reads` function supports a `filter_clean` parameter. When enabled, it strictly filters out any read records where the target junction is flanked by sequence variants. This includes both split-reads and mate-pair bridges.
