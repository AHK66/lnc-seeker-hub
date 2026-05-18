# lnc_seeker_bokeh Architecture

This document describes the modular architecture of the `lnc_seeker_bokeh` package, which provides a Bokeh-based interactive interface for analyzing genomic BAM data and BAM-derived caches processed by the Rust backend.

## Design Philosophy

The app follows an **Orchestrator-Manager pattern** to maintain a clean separation of concerns. It is designed for **BAM-free operation**, prioritizing high-speed hydration from `.lnc_cache.bin` files to allow visualization without access to original BAM files.

## Core Components

### 1. Orchestrator
*   **[ui_manager.py](ui_manager.py)**: The central `VisualizerApp` class. It manages the Bokeh document life-cycle, initializes all sub-managers, and coordinates the high-level update loop (`update_all_samples`).

### 2. UI & Layout
*   **[ui_layout.py](ui_layout.py)**: Handles the instantiation of Bokeh widgets and defines the initial structural layout. This decouples the visual definition from the behavioral logic.
*   **[constants.py](constants.py)**: Centralizes CSS styles, HTML templates for progress bars, and color palettes (e.g., `RedGrayBlue11`).

### 3. Feature Managers
These classes handle the logic for specific interactive sections of the application:
*   **[transcript_creator.py](transcript_creator.py)**: Manages the "Transcript Creator" tab, allowing users to manually assemble transcripts from detected junctions.
*   **[reads_manager.py](reads_manager.py)**: Orchestrates the detailed "Supporting Reads" view, including fetching read alignments for specific junctions from the **BAM-derived cache**.
*   **[genome_manager.py](genome_manager.py)**: Manages the genome annotation track (exons, introns, CDS).
*   **[rules_manager.py](rules_manager.py)**: Handles the dynamic generation of comparative highlighting rules based on sample groups.
*   **[selection_manager.py](selection_manager.py)**: Synchronizes junction selections across different plots to ensure a unified user experience.

### 4. Plotting & Rendering
*   **[coverage_plot.py](coverage_plot.py)**: Contains the lifecycle logic for sample-specific coverage tracks, including data updates and figure resets.
*   **[plotting_base.py](plotting_base.py)**: Low-level utility functions for visual styling, color mapping, and tool configuration.

## Technical Implementation Details

For in-depth technical specifications and logic descriptions, see [TECHNICAL_DETAILS.md](TECHNICAL_DETAILS.md). Key topics include:
*   **Automatic Coverage Resampling**: Logic for dynamic, cliff-aware depth visualization during zoom/pan.
*   **The Flat Model**: Logic for merging features into physical footprints and deriving introns.
*   **Session Isolation**: The "Mailbox" approach for thread-safe UI updates.
*   **Concurrency Control**: The "Busy Gate" mechanism for serialized backend execution.

## Data & State
*   **[shared_data.py](shared_data.py)**: Defines the session-specific `state` schema. Every browser tab/session receives a fresh, isolated state dictionary.
*   **[state.py](state.py)**: Contains stateless utility functions for logging, environment verification, and configuration loading. All functions require an explicit `state` object.
*   **[data_utils.py](data_utils.py)**: Contains heavy lifting for data transformation, such as parsing annotation dataframes, calculating global coordinate ranges, and determining "Marked Sets" for comparative analysis.
*   **[pipeline.py](pipeline.py)**: Orchestrates background analysis threads. It bridges Python and the Rust backend using the Mailbox approach to safely update session-specific UI.

## Data Flow

1.  **Input**: User selects a gene or samples. Selection is stored in the session-specific `state`.
2.  **Trigger**: `ui_manager` executes `run_analysis_thread(state)`.
3.  **Rust Backend**: The thread calls the stateless Rust library (`lnc_seeker`). Progress is tracked via an atomic `SessionProgress` object passed from the Python state.
4.  **Handoff**: Once Rust finishes, the background thread updates `state["analysis_data"]` and calls `doc.add_next_tick_callback`.
5.  **Update**: On the next event loop tick, the UI thread runs `refresh_session_ui`, triggering all managers to update their `ColumnDataSource` objects from the isolated `state`.
6.  **Rendering**: Bokeh synchronizes changes to the specific browser tab that initiated the request.
