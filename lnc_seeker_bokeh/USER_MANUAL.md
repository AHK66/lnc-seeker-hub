# lncRNA Seeker Hub - User Manual

The **lncRNA Seeker Hub** is a Bokeh-based interactive visualizer designed for high-performance genomic analysis. It utilizes a "BAM-free" architecture, prioritizing speed by hydrating data from pre-computed caches (`.lnc_cache.bin`) while allowing real-time interaction and deep-diving into read-level evidence.

---

## 📖 Quick Start

1. **Launch the App**: Run the Bokeh server as described in the [Quickstart](../README.md#quickstart).
2. **Select Gene**: Use the sidebar to choose a gene for analysis.
3. **Choose Samples**: Select one or more samples/cohorts from the `Data Selection` panel.
4. **Interact**: Use the mouse to zoom, pan, and hover over junctions for details.

---

## 🛠️ Main Panels

### 1. Analysis Controls (Sidebar)

The primary interface for data selection and global plot settings.

- **Gene & Sample Selection**: Select the target gene and samples.
- **Data Selection**: Toggle specific GTF files to be used as reference annotations.
- **Coverage & Junctions Profile**: Adjust track height, mapping-quality filtering, background vs. foreground display, and normalization behavior.
- **Junction Filtering**: Filter by minimum read support or mismatch/insertion flanks.

### 2. Supporting Reads (Full Read Layout)

Provides a detailed, alignment-level view of individual reads.

- **Read Visualization**: Inspect individual read alignments, insertions (`I`), and mismatches (`X`).
- **Gap & Deletion Visuals**: Toggle markers for deletions and numerical labels for intron gap sizes.

### 3. Transcript Creator

Allows users to manually assemble and export new transcript models.

- **Manual Assembly**: Add selected junctions from the main plot to the "Creator" list.
- **Refinement**: Manually set start/end positions, strand, and IDs.
- **Export**: Save models as **GTF**, **GFF3**, or **JSON** for later re-import.

## 🖱️ Interactive Features

### Unified Genomic View (Sync Crosshair)

- **Synchronized Zoom/Pan**: Zooming into a region on the coverage plot automatically updates the Annotation and Read plots (`Ctrl + Mouse Wheel`).
- **Sync Crosshair**: A vertical line that follows the mouse across all plots for coordinate-perfect inspection. Toggle in `General` settings.

### Comparative Highlighting

Define "Rules" to highlight junctions based on their presence in different sample sets (e.g., "Show junctions only present in Sample A but not Sample B").

---

## 📍 Coordinate System

The application uses a **zero-based, half-open** interval system $[start, end)$.

- **Example**: An exon from $100$ to $200$ (0-based) translates to $[101, 200]$ in 1-based inclusive tools like IGV.

---

## ⚙️ Configuration & Performance

- **Output Backends**: Toggle between `Canvas` (performance) and `SVG` (publication quality) in `General` settings.
- **Session Import**: Use the `File Import (JSON)` button in the Transcript Editor to restore previous sessions.
- **Config**: Start from [config.template.json](../config.template.json) and keep your local runtime settings in `config.json`.

For deeper technical details, refer to the [Technical Details](TECHNICAL_DETAILS.md) and [Architecture Overview](ARCHITECTURE.md).
