# lncRNA Seeker Hub

A small toolkit (Rust backend + Python UI) for interactive analysis and visualization of splice junctions and coverage from BAM-derived caches. The UI is a Bokeh application (lncRNA Seeker Hub) backed by a Rust analysis core compiled as a Python extension.

## Table of contents

- Overview
- Features
- Architecture
- Installation
- Quickstart
- BAM Collection & Cache Tool (lnc-seeker-collect)
- Configuration
- Detailed Documentation Index
- Development
- Contributing
- License & Authors
- Contact

## Overview

lnc-seeker performs coverage and splice-junction analysis on BAM files and provides an interactive Bokeh dashboard. It utilizes a **portable compact storage format** (`.lnc_cache.bin`), creating a **BAM-derived cache** that allows the UI to perform analysis and visualization even when the original BAM files are no longer available.

Target audience: bioinformaticians and researchers analyzing splice junctions and local coverage patterns.

## Features

- **Fast Rust-based Processing**: High-performance backend using `noodles` for BAM files and high-speed multi-block hydration for caches.
- **BAM-Free Interactive Dashboards**: Bokeh application for deep exploration of coverage, junction arcs, and full-read layouts without needing original BAM files.
- **Portable Compact Storage**: Generated caches are fully self-contained (including reference metadata), allowing the UI client to function independently of the source BAMs.
- **Genome-Verified Mismatch Discovery**: Reliable identification of sequence mismatches by comparing alignments against an indexed reference genome (FASTA). Sequence variations are highlighted in the read layout with labeled boxes (`X` for mismatches, `I` for insertions) and deletion markers.
- **Cross-Platform Compatibility**: Standardized BGZF-compressed BAM output and unified cache resolution across Linux and Windows.
- **Config-driven analysis**: Flexible selection of BAMs or pre-computed caches and visualization options.
- **Parallel Cache Hydration**: High-speed multi-block cache loading (`LNC1` format) that utilizes multi-core parallelism to decompress metadata and read names (Differential Huffman) concurrently.

## Architecture

- Rust analysis core: [lnc_seeker_lib/src/lib.rs](lnc_seeker_lib/src/lib.rs) compiled into a Python extension via maturin / pyo3.
- BAM Collection & Cache Tool: [lnc_seeker_collect](lnc_seeker_collect/) standalone CLI application for subsetting BAM files and generating caches.
- Python UI: [lnc_seeker_server.py](lnc_seeker_server.py)  Bokeh server app.
- Build metadata: [pyproject.toml](pyproject.toml).

## Installation

### Prerequisites

- **Python 3.8+**: Ensure you have the ability to create virtual environments. On Debian/Ubuntu systems, you might need to install the `python3-venv` package (e.g., `sudo apt install python3-venv`).
- **Rust Toolchain**: Required to build the high-performance analysis core. [Install rustup](https://rustup.rs/).

### Setup Steps

1. **Clone the repository:**

   ```powershell
   git clone https://github.com/AHK66/lnc-seeker-hub.git
   cd lnc-seeker-hub
   ```

2. **Create and activate a virtual environment:**

   ```powershell
   # Windows (PowerShell)
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1

   # Unix/macOS
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Install dependencies:**

   ```powershell
   # Windows
   python -m pip install -r requirements.txt
   python -m pip install -r requirements-build.txt

   # Unix/macOS
   python3 -m pip install -r requirements.txt
   python3 -m pip install -r requirements-build.txt
   ```

4. **Build the Rust extension:**
   This project uses `maturin` to compile the Rust core as a Python module (`lnc_seeker`).

   ```powershell
   # Windows
   python -m maturin develop --release

   # Unix/macOS
   python3 -m maturin develop --release
   ```

## Quickstart

Once the installation is complete, you can launch the interactive dashboard:

1. **Configure**: Copy [config.template.json](config.template.json) to `config.json`, then list your BAM files or paths to existing `.lnc_cache.bin` caches in the local `config.json`.
2. **Launch**:

   ```powershell
   python -m bokeh serve --show lnc_seeker_server.py --session-token-expiration 3600
   ```

## BAM Collection & Cache Tool (lnc-seeker-collect)

The toolkit includes a standalone Rust application for subsetting large BAM files and generating the optimized `.lnc_cache.bin` format used for visualization. It extracts records matching specific gene regions or custom coordinates.

### Build

```powershell
cargo build --release --bin lnc-seeker-collect
```

### Usage

```powershell
./target/release/lnc-seeker-collect path/to/config.cfg
```

A sample configuration can be found in [lnc_seeker_collect/config.cfg](lnc_seeker_collect/config.cfg).

## Configuration

Primary local config file: `config.json` (kept out of version control).
Use [config.template.json](config.template.json) as the canonical public example, then copy it locally to `config.json` and fill in your dataset-specific paths.

Key sections:

- `data_selection`: select data sources via `bam_paths` (supports BAM files, `.lnc_cache.bin` paths, or a mapping gene -> [paths]), `filter_outliers`, and an optional `genome_path` (indexed FASTA) for mismatch verification. For **BAM-free operation**, provide the paths to your pre-computed caches here.
- `coverage_and_junctions_profile`: plotting heights, mapping-quality thresholds, coverage style options.
- `junctions_splicing`, `genome_annotations`, `full_read_layout`, `general`: UI visual and behavior options.

See [config.template.json](config.template.json) for the full public example.

## Detailed Documentation Index

Detailed technical documentation for the various components of the toolkit:

- **General & Infrastructure**
  - [Memory Architecture](MEMORY_ARCHITECTURE.md): Details on telemetry metrics, memory management strategies, and optimization patterns for the Rust-Python bridge.
  - [Apache Configuration](APACHE_CONFIG.md): Technical guide for configuring an Apache 2 server as a reverse proxy for the Bokeh application.
- **Interactive Dashboard (lnc_seeker_bokeh)**
  - [User Manual](lnc_seeker_bokeh/USER_MANUAL.md): High-level guide to the interactive dashboard, panels, and interdisciplinary controls.
  - [Architecture Overview](lnc_seeker_bokeh/ARCHITECTURE.md): Modular "Orchestrator-Manager" pattern, session management, and UI logic.
  - [Technical Details](lnc_seeker_bokeh/TECHNICAL_DETAILS.md): Specifications for features like the "Cliff-Aware Resampling Algorithm" and concurrency control.
- **Data Collection (lnc-seeker-collect)**
  - [Extraction Process](lnc_seeker_collect/PROCESS.md): Algorithmic details on read remapping, stable region hashing, and extraction logic.
  - [Configuration Format](lnc_seeker_collect/CONFIG_FORMAT.md): Detailed syntax guide for the `lnc-seeker-collect` configuration tool.
- **Core Library (lnc-seeker-lib)**
  - [Library Overview](lnc_seeker_lib/src/OVERVIEW.md): High-level orientation of the Rust core architecture and PyO3 C-API.
  - [Annotation Logic](lnc_seeker_lib/ANNOTATION_LOGIC.md): Multi-pass GTF retrieval and coordinate normalization for full transcript integrity.
  - [Cache Architecture](lnc_seeker_lib/CACHE_ARCHITECTURE.md): Specification of the portable compact storage format (`.lnc_cache.bin`).
  - [Header Compression](lnc_seeker_lib/HEADER_COMPRESSION.md): Details on the dictionary-based Huffman scheme for BAM read names.
  - [Test Documentation](lnc_seeker_lib/TESTS.md): Overview of unit and integration tests within the library.

## Development

- Build: `python -m maturin develop` (see [Installation](#installation)).
- Runtime Python dependencies: `requirements.txt`.
- Build tooling for the Rust extension: `requirements-build.txt`.
- Formatting & linting: follow Python and Rust conventions (black / rustfmt).
- Tests: Rust unit and integration coverage lives in [lnc_seeker_lib/TESTS.md](lnc_seeker_lib/TESTS.md); run `cargo test` from `lnc_seeker_lib` for the current fast suite. Python-side pytest coverage is not added yet.

Notes:

- The Rust API accepts both an array and a mapping form for `data_selection.bam_paths` (the mapping is flattened internally).
- If Pylance/VS Code cannot resolve the compiled module, select the workspace `.venv` interpreter and rebuild the extension with `python -m maturin develop`.

## Contributing

- Fork, branch, and open a PR with focused changes.
- Include tests where appropriate and update documentation for new features.

## License & Authors

**Authors:** Arne Kutzner and Pok-Son Kim

This project is licensed under the Apache License, Version 2.0. See the [LICENSE](LICENSE) file for details.

## Contact / Support

Open issues on GitHub or contact the maintainers listed in the repository metadata.
