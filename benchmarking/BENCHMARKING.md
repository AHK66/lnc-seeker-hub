# Benchmarking & Performance Tracking

This document explains how to use the built-in benchmarking toolkit in **lnc-seeker-hub** to measure performance, memory consumption, and processing latency.

## Overview

The benchmarking system consists of three main components:
1.  **`BenchmarkManager` (Python)**: A singleton manager in [benchmarking/python/benchmark_manager.py](benchmarking/python/benchmark_manager.py) that logs events to a CSV file.
2.  **`BenchmarkMonitor` (Rust)**: High-precision memory and duration tracking in the Rust core ([benchmarking/rust/benchmark_core.rs](benchmarking/rust/benchmark_core.rs)).
3.  **Visualization Script**: A Bokeh-based reporter that generates interactive HTML dashboards from the logs.

## 1. Enabling Benchmarking

Benchmarking is disabled by default to minimize performance overhead. To enable it:

1. Copy [config.template.json](../config.template.json) to `config.json` if you have not created a local config yet.
2. Open the local `config.json`.
3. Locate the `general` section.
4. Set `enable_benchmarking` to `true`:

```json
"general": {
    "enable_benchmarking": true
}
```

When enabled, the application will record performance metrics for various operations (initialization, data processing, UI interactions) and store them in [benchmarking/data/performance_log.csv](benchmarking/data/performance_log.csv).

## 2. Running a Benchmarking Session

Once benchmarking is enabled in the local configuration:

1. Start the Bokeh server as usual:
   ```powershell
   python -m bokeh serve --show lnc_seeker_server.py
   ```
2. Interact with the dashboard (e.g., select genes, change coverage settings, toggle mismatch labels).
3. The logs are updated in real-time. Each event records:
    - **Duration**: Time spent on the operation.
    - **Memory (RSS)**: Total process memory usage.
    - **Cache Details**: Breakdown of memory used by the Rust core, annotations, and JSON serialization.

## 3. Generating the Performance Report

After collecting data, you can generate an interactive HTML report:

```powershell
python benchmarking/scripts/visualize_performance.py
```

This script processes the CSV logs and creates a self-contained report at [benchmarking/data/performance_report.html](benchmarking/data/performance_report.html).

### Report Contents
- **Activity Timeline**: A Gantt chart showing the sequence and duration of processing events.
- **Memory Profile**: A multi-layered area chart showing how different components (Core, Annotations, JSON Overhead) contribute to the total memory footprint.
- **Latency Analysis**: Statistics on time spent in various stages of the pipeline.

## 4. Technical Details

### Python Implementation
The `BenchmarkManager` uses `psutil` for system-level memory monitoring and provides decorators to wrap UI callbacks and processing pipelines.

### Rust Implementation
The `BenchmarkMonitor` provides atomic peak memory tracking and high-resolution timing, which is bridged to Python via PyO3.

## Supporting Files
- **Log Data**: [benchmarking/data/performance_log.csv](benchmarking/data/performance_log.csv)
- **Report**: [benchmarking/data/performance_report.html](benchmarking/data/performance_report.html)
- **Label Config**: [benchmarking/scripts/visualize_performance_labels.json](benchmarking/scripts/visualize_performance_labels.json) (Customizes report labels)
