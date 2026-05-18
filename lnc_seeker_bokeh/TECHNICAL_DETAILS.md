# lnc_seeker_bokeh Technical Notes

This document contains detailed technical specifications and logic descriptions for the `lnc_seeker_bokeh` package.

## Automatic Coverage Profile Computation

The coverage profile is a sample-specific visualization that dynamically computes and renders depth data from BAM files or pre-computed caches, adapting to the user's current zoom level to maintain high performance without sacrificing visual accuracy.

### 1. Data Source
Raw depth data (including all reads and high-quality reads) is calculated by the Rust backend. It is stored as NumPy arrays in the session-specific `state["analysis_data"]`. The system tracks two distinct profiles:
- **Background (`y_bg`)**: Represents all reads.
- **Foreground (`y_fg`)**: Represents filtered high-quality reads.

**BAM-Free Operation**: When loading data, the backend prioritizes `.lnc_cache.bin` files (the **BAM-derived cache**). If these are provided, the UI can render all coverage and junction data without the original BAM files being present, enabling fast, portable analysis.

### 2. Dynamic Range Adaptation (Zoom Logic)
The app listens for changes to the shared genomic range (`shared_x_range.start` and `shared_x_range.end`). To ensure UI responsiveness during smooth panning or zooming, the resampling is **debounced by 200ms**. When the range settles, the `on_range_change` callback triggers `update_sample_coverage` in [lnc_seeker_bokeh/coverage_plot.py](coverage_plot.py) for all active sample plots.

### 3. Cliff-Aware Resampling Algorithm
The core "automatic coverage profile computation" follows these steps:
1.  **Resolution Detection**: It determines the `pixel_width` of the target diagram (typically the `inner_width` of the Bokeh figure, defaulting to 1200 pixels).
2.  **Viewport Filtering**: Slices the global depth arrays to the current visible genomic domain, plus a small margin (50 points) to prevent artifacts at the edges.
3.  **Resampling Logic**:
    - If the number of genomic points in view is less than or equal to the `pixel_width`, the raw data is used.
    - If it exceeds the `pixel_width`, it calls a Rust binding (`lnc_seeker.downsample_coverage_py`).
4.  **Preserving Junctions (Cliffs)**: The resampler is "cliff-aware." It receives the positions of all detected splice junctions (cliffs) within the current view. The algorithm ensures that these critical transition points and depth peaks are preserved in the downsampled output, preventing the "smoothing out" of biological signals that often happens with simple averaging or decimation.
5.  **UI Update**: The downsampled data is pushed to the `ColumnDataSource`, triggering an efficient update of the `step` glyphs in the browser.

## Annotation Processing & Intron Detection

The app transforms raw GTF/GFF records into a structural reference used to categorize splice junctions detected in BAM files.

### 1. The Flat Model (Physical Footprints)
To accurately identify introns, the lnc_seeker_bokeh builds a "Flat Model" in [lnc_seeker_bokeh/data_utils.py](data_utils.py). 
*   **Merging Strategy**: For each transcript ID, all structural features (`exon`, `CDS`, `UTR`, `ncRNA`) are merged into contiguous physical blocks. 
*   **Intron Derivation**: True biological introns are defined as the **gaps** between these merged blocks.
*   **Redundancy Handling**: Merging is critical because GTFs often contain redundant or overlapping records (e.g., an `exon` record and a `CDS` record for the same range). Without merging, the gap-detection logic would produce "zero-length introns" and fail to match BAM junctions.

### 2. Feature Filtering
The selection of features for the Flat Model is intentionally restrictive:
*   **Excluded Summary Features**: Records labeled as `transcript` or `gene` are **filtered out** before merging. Because these records span the entire length of the gene, including them would cause the merging algorithm to "swallow" all introns, merging the entire gene into a single giant block.
*   **Included Structural Features**: Only features that define physical boundaries (like `exon`, `CDS`, and `UTR`) are used to calculate the backbone and intron sets.

### 3. Junction Categorization
Splice junctions detected in BAM files or retrieved from sub-sampled caches are compared against the Flat Model of the reference genome:
*   **Curated**: The junction's start and end coordinates exactly match an intron gap from a curated transcript (e.g., `NM_` or `NR_` prefixes).
*   **Predicted**: The junction matches an intron from an uncurated or computational transcript (e.g., `XM_` or `XR_`).
*   **Novel**: The junction does not exist in the reference Flat Model.

### 4. Coordinate Consistency
Genomic data processing requires strict adherence to coordinate systems:
*   **0-based vs 1-based**: The system normalizes GTF (1-based) coordinates to match BAM (0-based) alignment starts.
*   **Offsets**: When a `gtf_offset` is provided in the configuration, it is applied **only to the start coordinate** of every feature. This is intentional behavior designed to reconcile specific GTF coordinate conventions with the 0-based system used by the BAM backend.

## Junction & Read Filtering (Mismatches/Insertions)

To improve the detection of high-confidence splice junctions, `lnc_seeker` implements a filtering system that identifies and excludes junctions flanked by alignment artifacts or biological variants (mismatches and insertions).

### 1. Algorithmic Principle: Flank Quality Control
Splice junctions are defined by the transition between two alignment segments separated by a "Skip" (`N`) CIGAR operation. The system evaluates the **quality of the anchor points** immediately surrounding these transitions:
- **Clean Junction**: A junction where both the upstream segment (ending at the splice site) and the downstream segment (starting after the splice site) have **zero mismatches** and **zero insertions**.
- **Excluded Junction**: A junction where either flanking segment contains at least one mismatch or insertion.

### 2. Implementation Levels

#### A. Backend Aggregation (`coverage.rs`)
During the initial BAM analysis, the `JunctionStore` tracks two counters for every unique coordinate pair:
- `total_reads`: All reads supporting the junction.
- `reads_clean`: Only those reads passing the Clean Junction criteria.

#### B. UI Counter Labels (`coverage_plot.py`)
When the "Filter Mismatch/Insertion Flanks" option is toggled, the coverage plot's junction arcs dynamically switch their labels. The labels are updated to show only the `reads_clean` count, and the arc thicknesses are scaled accordingly.

#### C. Full Read Filtering (`processing.rs`)
The filtering is most aggressive in the "Full Reads Layout" diagram. When activated:
1. **Target Filtering**: The backend re-validates supporting reads for the selected junction. If a read supports the target junction but that specific junction instance is "dirty" in that read, the **entire read** is removed from the dataset.
2. **Internal Junction Hiding**: If a read is clean at the *target* junction and thus visible, the rendering process in [lnc_seeker_bokeh/reads_manager.py](reads_manager.py) further checks all other internal introns. Any internal junction that is "dirty" will have its intron line hidden, leaving only the clean backbone visible.

### 3. Bridged Junctions (Mate-Pair Support)
For junctions formed by the physical gap between two mates (where the junction itself is not directly sequenced in a single read but implied by the pairing), the system applies the same strictness:
- The segment of **Mate 1** that terminates at the junction boundary must be clean.
- The segment of **Mate 2** that initiates at the junction boundary must be clean.
- If either mate's anchor segment contains variants, the bridged junction is considered "dirty."

## Advanced Concurrency & Race Condition Mitigation (The "Coalescing Gate")

To ensure robust behavior during rapid user interactions (e.g., quickly switching between genes or changing settings while a scan is active), the pipeline employs a sophisticated synchronization layer:

1.  **Atomic State Transitions (Finalization Order)**: In the [pipeline.py](pipeline.py) orchestrator, the `analysis_running` flag is transitioned to `False` **before** the UI notification callbacks are executed. This prevents race conditions where the UI thread might poll for status and see a "running=True" state even though the data is technically ready, which previously caused diagrams to occasionally "vanish" or fail to render the final iteration.

2.  **Request Coalescing**: Instead of simply ignoring new requests while the Rust backend is busy, the orchestrator now employs a coalescing loop. If parameters change while a scan is underway, the current thread will detect the configuration mismatch upon completion and **immediately restart** the analysis. This ensures that the UI always reflects the *latest* user intent without spawning multiple competing threads.

3.  **Staleness Guard (Post-Rust Validation)**: Heavy-duty Rust analysis can take several seconds. To prevent stale data from being "injected" into a newer session state (e.g., if a user switched genes while Rust was busy), the pipeline performs a validation check immediately after the Rust binding returns. If the configuration has been superseded, the results are discarded before they can reach the state's `analysis_data`.

4.  **Cross-Gene Integrity Control (`data_gene_name`)**: The system tracks the `data_gene_name` associated with the current `analysis_data`. The UI manager verifies this mapping before rendering any experimental content. This serves as a final shield against "phantom diagrams" where plots from a previous gene selection might briefly appear in the layout of a newly selected gene due to asynchronous update delays.

5.  **Magic-Key (Request IDs) for Fetches**: The [ReadsManager](reads_manager.py) uses an incrementing `request_id` to version asynchronous BAM fetch operations. Only the callback associated with the *most recent* request is permitted to update the Bokeh `ColumnDataSource` or clear the progress status, solving the problem of persistent "Fetching reads..." messages caused by out-of-order network/disk returns.

## Session Isolation (The Mailbox Approach)

To ensure stability in multi-user environments, the lnc_seeker_bokeh uses an **Explicit State / Mailbox** architecture:

1.  **Isolation**: There are **zero globals**. Every Bokeh session (`VisualizerApp`) owns a private `state` dictionary initialized via `shared_data.create_session_state()`.
2.  **Explicit Passing**: The `state` object is passed into every function and sub-manager. Components do not use `curdoc()` or global registries to find data.
3.  **Background Processing**: Heavy tasks (like the Rust analysis or depth calculation) run in standard Python `threading.Thread` instances.
4.  **The Mailbox (Synchronization)**: Since background threads cannot safely modify Bokeh models directly, they "post" updates back to the UI thread using `doc.add_next_tick_callback(callback)`. The UI thread then picks up these "messages" from the session's state and applies them to the renderers.
