# Memory Architecture and Benchmarking Report

This document details the memory management strategies, telemetry metrics, and optimization patterns used in the `lnc-seeker` project to bridge the high-performance Rust backend with the interactive Python/Bokeh frontend.

## 1. Benchmarked Memory Categories

The performance telemetry suite monitors the Resident Set Size (RSS) of the application and decomposes it into specific architectural components.

| Component | Logic | Description |
| :--- | :--- | :--- |
| **Initial OS RSS** | `baseline` | The memory footprint of the Python interpreter and Bokeh server upon initialization. |
| **Cache Core (LNC1)** | `Rust-owned` | Memory allocated within the Rust backend for indexing and caching high-depth sequencing data (BAM/CRAM). |
| **Cache Annotation** | `Rust-owned` | Memory allocated within Rust for the optimized LBA (Low-latency Binary Annotation) cache and GTF indices. |
| **JSON/PyObject Bridge** | `Transient` | The raw size of the UTF-8 JSON string used to serialize data between Rust and Python. |
| **Cache Related** | `Mixed` | Transient memory for supplementary data like active coverage vectors or junction points currently being processed. |
| **Peak Processing** | `Python-owned` | **The "Black Area"**: Calculated as `Total RSS - Measured Components`. It represents the overhead of Python objects, Pandas DataFrames, and internal VM allocations. |

---

## 2. The JSON Overhead and Serialization Bridge

### The "JSON Explosion" Problem

`lnc-seeker` uses a JSON bridge to pass complex analysis results from Rust to Python. While convenient, this creates a significant memory bottleneck:

1. **String Duplication**: Rust generates a massive string; Python holds that string in memory.
2. **Object Inflation**: `json.loads()` converts the string into millions of individual Python objects (`dict`, `list`, `float`, `str`). Each annotation feature can inflate from ~40 bytes in Rust to ~400 bytes in Python.
3. **Immutability**: Python strings are immutable; any intermediate manipulation creates new copies.

### Alternative Approaches & Disadvantages

| Approach | Description | Disadvantages |
| :--- | :--- | :--- |
| **Apache Arrow** | Shared memory buffers with zero-copy access. | Extremely high implementation complexity; requires strict memory alignment and schema management. |
| **Protocol Buffers** | Binary serialization format. | Requires a complex compilation step; while faster than JSON, it still incurs deserialization overhead into Python objects. |
| **Raw C-FFI Buffers** | Passing raw pointers to Rust-allocated arrays. | Highly unsafe; risks segmentation faults if the Rust lifecycle isn't perfectly synchronized with Python's Garbage Collector. |
| **LBA Spatial Queries** | Rust only returns features within the current viewport. | **The "Disadvantage of Interactivity"**: Panning the plot becomes jittery as every move requires a new backend request instead of a smooth local UI filter. |

---

## 3. Python VM and Garbage Collection Optimization

### The "Black Plateau" Phenomenon

In high-throughput biological tools, the Python Virtual Machine often exhibits a "Black Plateau" in memory diagrams. This is caused by the Python Garbage Collector (GC):

* **Deferred Reclamation**: Python uses reference counting and a cyclic GC. Even after a large object (like the raw JSON string) is no longer used, the GC may not immediately return those memory pages to the Operating System.
* **Fragmentation**: Allocating millions of small objects (annotations) fragments the heap, making it difficult for the OS to reclaim space.

### Mitigation: The "Forced Early Drop" Strategy

To stabilize memory during long sessions, `lnc-seeker` implements a "Forced Early Drop" pattern in the [data pipeline](lnc_seeker_bokeh/pipeline.py):

1. **Redundancy Pruning**:
   As soon as a JSON payload is converted into a Pandas DataFrame, the project explicitly prunes the original dictionary:

   ```python
   state["analysis_data"].pop("annotations", None) # Remove raw list
   ```

2. **Explicit Deletion**:
   Using `del` on the raw string and intermediate dictionaries to drop reference counts to zero immediately.

   ```python
   del res_json # Destroy the raw transfer string
   ```

3. **Manual Collection**:
   Calling `gc.collect()` at the end of a processing cycle. While usually discouraged in general Python scripts, it is essential in interactive apps to Recede the "Black Plateau" and provide a stable memory ceiling for the user.
4. **Linux Heap Trimming**:
   On Linux, the pipeline follows `gc.collect()` with `malloc_trim(0)` so glibc can return free heap pages to the OS more aggressively after large deserialization and Pandas-heavy processing bursts.

### Disadvantage of Forced Dropping

* **CPU Latency**: Calling `gc.collect()` is a "stop-the-world" event for the Python interpreter. It incurs a minor latency penalty (~50-200ms) at the end of every data update, trading a small amount of UI smoothness for significant memory reliability.
* **Platform Variance**: `malloc_trim(0)` is Linux/glibc-specific, so it can improve RSS reclamation on Linux while having no effect on Windows and little effect when fragmentation is dominated by long-lived objects.
