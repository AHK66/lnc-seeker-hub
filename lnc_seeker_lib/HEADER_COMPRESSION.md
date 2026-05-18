# Experimental Header Compression (Differential Huffman)

This document describes the experimental compression scheme used for BAM QNAMEs (read names) in the `lnc_seeker_hub` cache system.

## Motivation
BAM read names often exhibit high redundancy, especially when sorted or localized to specific genomic regions (e.g., Illumina flowcell coordinates or Nanopore run IDs). Traditional Zstd compression on the whole cache is effective but can be slow to decompress during interactive visualization. This scheme provides a domain-specific compression that prioritizes **decompression speed** while achieving significant space savings.

## 1. Differential Encoding (Delta)
The core of the compression is a delta-over-predecessor approach. Instead of storing each string, we store the edits required to transform the previous name into the current one.

### Edit Operations (`EditOp`)
We use a prefix-biased diffing algorithm that generates four types of operations:
- **Match(N)**: Keep N characters from the previous string (starting from the beginning).
- **Delete(N)**: Skip N characters of the previous string.
- **Insert(String)** / **Substitute(String)**: Add or replace bytes with a literal string.

By focusing on common prefixes, we exploit the Illumina "coordinate" structure where many reads share the same lane/tile info.

## 2. Huffman Coding
The stream of edit operations and literal bytes is entropy-coded using **Canonical Huffman Coding**.

### Symbol Space
The encoder uses a dynamic 16-bit symbol space to balance efficiency and complexity:

| Symbols | Operation | Range | Details |
|---------|-----------|-------|---------|
| 0-255 | Literal | Single Byte | Raw ASCII value |
| 256-511 | Match | 1-255 | Bytes to copy from predecessor |
| 512-767 | Delete | 1-255 | Bytes to skip in predecessor (offset) |
| 768-1023 | Substitute | 1-255 | **Optional**: Bytes to replace (advances offset) |
| 768 / 1024 | End Symbol | N/A | Marks the end of a single read name |

#### Two-Mode Alphabet
To optimize metadata overhead, the library supports two modes:
1. **Compact Mode** (Default): Disables `Substitute` ops. The alphabet size is **769 symbols** (End symbol at 768). This is faster and uses slightly less header space.
2. **Extended Mode**: Enables `Substitute` ops. The alphabet size is **1025 symbols** (End symbol at 1024). This provides better compression for names with small internal variations.

This setting is controlled via `header_compression_use_substitutes` in `config.cfg`.

### Canonical Format
The compression uses the Canonical Huffman property, where only the list of code lengths is required to reconstruct the decoding table. This minimizes metadata overhead.

## 3. Data Structure
The compressed headers are stored in the `CompressedHeaders` struct:
- `lengths: Vec<u8>`: The 770 code lengths for the Huffman tree.
- `bitstream: Vec<u8>`: The actual compressed data.

## 4. Performance & The "Zstd Architectural Trap"
- **The Trap**: Compressing an entropy-coded bitstream (like Huffman data) inside a general-purpose Zstd stream creates a serial dependency. The CPU must wait for Zstd to finish before Huffman can begin.
- **Decompression**: Optimized for $O(1)$ table-driven decoding. It is designed to be significantly faster than general-purpose compression when hydrating tens of thousands of read names for the UI.
- **Parallelism**: By storing headers in a dedicated block outside the main Zstd payload, we can utilize `rayon::join` to decompress metadata and read names concurrently, effectively hiding the decompression latency.
- **Single-Pass**: The decoder uses a high-performance single-pass implementation that avoids intermediate vector allocations, maintaining hydration times under 1 second for typical datasets.
- **Memory**: Once decompressed in RAM, the strings are stored in a standard `Vec<String>`. The library ensures that the compressed bitstream is cleared from memory after hydration to avoid double-counting.

## 5. Integration
The compression is triggered via the `--compress-headers` flag in `lnc-seeker-collect`. It populates the `compressed_headers` field in `BamCache`, which is automatically detected and decompressed by the Rust `ReadProvider` whenever a cache file is loaded.
