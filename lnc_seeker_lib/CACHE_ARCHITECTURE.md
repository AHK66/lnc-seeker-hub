# Cache Architecture: Tiered BAM-Derived Cache System

The `lnc_seeker_lib` implements a sophisticated, multi-tier caching layer (the **BAM-derived cache**) designed to support millions of genomic reads within a limited RAM footprint. This document describes the storage strategy, memory management, and persistent retrieval logic that enables **BAM-free workflows**.

## 1. Compact Flat Storage

To avoid the massive heap overhead of millions of small objects (like `String` read names or individual `Vec` segments), the cache uses a pooled storage architecture.

### Data Structures
- **`CompactRead`**: A 24-byte struct representing a single alignment.
    - Instead of a `String` name, it uses a `u32` name_id.
    - Instead of a `Vec<Segment>`, it uses a `u32` seg_offset and `u8` seg_count into a global segment pool.
- **`CompactSegment`**: A simple 8-byte `(start, end)` pair.
- **`segment_tags`**: A bit-packed `u8` per segment:
    - **Bit 0**: `is_followed_by_deletion` (boolean).
    - **Bits 1-4**: `mismatches` (count 0-15). If a reference genome is provided, these are verified against the sequence.
    - **Bits 5-7**: `insertions` (count 0-7).
- **`BamCache`**: The primary container for a single **BAM-derived cache** in RAM.

```rust
pub struct BamCache {
    pub reads: Vec<CompactRead>,      // Flat alignment storage
    pub segments: Vec<CompactSegment>, // Shared segment pool for this cache
    pub segment_tags: Vec<u8>,         // Bit-packed segment metadata
    pub names: Vec<String>,            // Unique name pool
    pub ref_names: Vec<String>,        // Reference sequence names
    pub ref_lengths: Vec<usize>,       // Reference sequence lengths (for BAM-free analysis)
    pub status: CacheStatus,           // Complete, Caching, or ExceededLimit
    pub compressed_headers: Option<CompressedHeaders>, // On-disk differential Huffman
    pub size_bytes: usize,             // Tracked memory footprint
}
```

## 2. Persistence & High-Speed Hydration (BAM-Free Workflows)

The persistence layer uses a custom multi-block binary format designed for maximum throughput, enabling **BAM-free workflows** where visualization requires only the pre-computed cache.

### The `LNC1` File Format
Caches are saved to `.lnc_cache.bin` files (the **BAM-derived cache**) using a multi-block structure:

1.  **Magic Header**: `LNC1` (4 bytes).
2.  **Compressed Header Block**: Differential Huffman-encoded read names. This block is stored outside the main Zstd stream to allow for specialized parallel decompression using `rayon::join`.
3.  **Metadata Block**: Zstd-compressed `Bincode` serialization of the remaining `BamCache` fields (read records, segments, indices).

### Parallel Loading Strategy
To minimize total load time, the `ReadProvider` implements parallel hydration using `rayon`:
- **Zstd Pipeline**: One task handles the decompression and deserialization of the large metadata block.
- **Huffman Pipeline**: A concurrent task runs the table-driven Huffman decoder to hydrate the `Vec<String>` from the bitstream.
- **Scaling**: This approach effectively hides the ~$1$s Huffman decompression time within the ~$6-7$s Zstd load time for large files (>1MB).

## 3. Differential Huffman Compression

Genomic read names have high redundancy (shared tile/lane info). A domain-specific compression scheme is used to minimize disk usage while outperforming general-purpose decoders:

1.  **Prefix-Biased Delta Encoding**: Only the differences between consecutive read names are stored (Match, Delete, Insert ops).
2.  **Canonical Huffman Coding**: The edit operations and literal characters are entropy-coded using $O(1)$ table-driven LUT decoding.
3.  **RAM Management**: Once hydratated into a `Vec<String>` in RAM, the compressed bitstream is automatically cleared to prevent double memory consumption.

## 4. Memory-Driven Lifecycle

Unlike traditional caches that use counts, this cache uses a **global memory budget** (configured via `max_cache_memory_mb`).

### Monitoring & Estimation
- **Incremental Estimation**: During the parallel analysis loop, threads track their estimated memory footprint incrementally (count based + object overhead). If the global sum of existing caches + the current estimate exceeds the limit, the operation is immediately aborted to prevent system instability.
- **Selection-Aware Purging**: When a new analysis starts, the cache manager purges all RAM records for samples that are not part of the current selection, keeping only the active cohort in memory plus any global background caches. Unselected samples remain available in their compressed `.lnc_cache.bin` form on disk.
- **Conservative Accurate Calculation**: When committing a cache, the system calculates `size_bytes` based on the `capacity()` of vectors and heap-allocated string buffers. 

## 3. Tiered Retrieval Logic

The architecture implements a streamlined retrieval strategy where analysis is decoupled from parsing, facilitating **portable, BAM-free workflows**:

1.  **Cache Population (Disk to RAM)**: 
    - If a sample is not in RAM, the system first checks for a `.lnc_cache.bin` file. 
    - **BAM-Free Logic**: If the input path ends in `.lnc_cache.bin`, the system loads it directly and performs analysis using internal metadata (reference names and lengths), completely bypassing any BAM file access.
    - If no cache exists and a `.bam` path is provided, it is populated via a fresh BAM scan.
2.  **Compact Analysis (RAM)**: All analysis (coverage, junctions) is performed exclusively from the `CompactRead` records and reference metadata in RAM. This ensures that any UI-driven recalculations (like MQ threshold changes) are virtually instant (~10-50ms per sample).
3.  **Persistence**: Caches are automatically serialized to `.lnc_cache.bin` once a sample analysis reaches the `Complete` status.

## 4. Path & Context Sensitivity

- **Normalization & Resolution**: All cache keys are normalized (lowercase and forward-slashes) to ensure consistency across Windows and Unix. The system automatically detects existing `.lnc_cache.bin` extensions to prevent redundant path decoration.
- **MQ Sensitivity**: Coverage calculations are performed on-the-fly. If the user changes the "Min Mapping Quality" in the UI, the system recalculates results from the `CompactRead` records in RAM, ensuring zero storage overhead for different filter configurations.

## 5. Thread Safety

The cache resides in a global `ReadProvider` managed by `OnceLock` and protected by an `RwLock<HashMap<...>>`. This configuration allows for:
- **Parallel Population**: Multiple BAM files being analyzed in parallel can write to the cache simultaneously.
- **Fast UI Retrieval**: The Python-bound Bokeh server can query the cache for filtered reads without blocking the analysis threads.

## 6. Persistent Disk Buffering

To ensure "instant-on" performance after application restarts, the library implements a background persistence layer.

### Binary Serialization (Bincode)
Instead of high-level formats like JSON, the cache uses **Bincode**. This is a compact binary format that maps directly to Rust struct memory layouts, eliminating the overhead of string parsing and number formatting.

### High-Speed Compression (Zstd)
Serialized data is compressed using **Zstd**. Zstd was chosen specifically for its ultra-fast decompression speeds, ensuring that the bottleneck for loading large samples is the disk I/O bandwidth rather than CPU cycles.

### Consistency & Safety
- **Memory Re-estimation**: When a binary cache is loaded from disk, the system re-calculates the local allocation overhead to ensure the global RAM budget tracking remains accurate.
- **Automatic Persistence**: Caches are automatically serialized to `.lnc_cache.bin` once a sample analysis reaches the `Complete` status.
