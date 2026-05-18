# TESTS for lnc_seeker_lib

This document describes the purpose and the properties asserted by the unit and integration tests in `lnc_seeker_lib`.

## Compression tests (src/compression.rs)

- **`test_all_compression_modes`**: Ensures header-name compression and decompression are a lossless round-trip for all supported `CompressionMode` values (`None`, `Huffman`, `Zstd`).
  - Purpose: validate `compress_header_block` and `decompress_header_block` behavior across modes.
  - Key assertions:
    - `compressed.num_names == names.len()` and `compressed.mode` matches the requested mode.
    - After decompression the returned vector of names equals the original input vector (exact string equality).
  - Properties tested: canonical Huffman encoding/decoding, zstd compression/decompression, and the raw (no compression) encoding path.

- **`test_empty_block`**: Verifies handling of empty input name lists.
  - Purpose: ensure `compress_header_block` and `decompress_header_block` handle zero-length inputs safely.
  - Key assertions:
    - `compressed.num_names == 0` and decompressed list is empty.
  - Properties tested: no-panics, correct length metadata, and empty-result semantics for Huffman mode.

## Cache & reads tests (src/tests.rs)

- **`test_cache_integrity`**: A unit test that exercises the in-memory/disk caching layer used by the reads provider.
  - Purpose: validate cache commit and subsequent retrieval, including reconstruction of split/paired records.
  - Setup: constructs mock `ReadInfo` records (one spliced read and one bridged pair) and calls `provider.commit_caching_compact(...)`.
  - Key assertions:
    - `get_filtered_reads(...)` returns `Some((filtered, status))` with `status == CacheStatus::Complete`.
    - The returned `filtered` set size matches expected reconstructed records (3 in the test: one spliced read + two records reconstructed from the bridged pair).
    - Returned record names include the expected read names (e.g. `read_A`, `read_B`).
  - Properties tested: cache write/commit, retrieval filtering by region, paired-read reconstruction logic, and cache status reporting.

- **`test_real_bam_cache`** (ignored by default): An integration test that runs against a real BAM file and exercises the full pipeline.
  - Purpose: compare cached retrieval performance and correctness against direct disk reads for a real junction.
  - Behavior:
    - Scans a specified BAM for the first spliced read to select a target junction.
    - Runs `run_analysis(...)` to populate the cache for that BAM.
    - Compares `get_junction_reads(...)` results from cache vs. disk and measures timings.
  - Key assertions:
    - Number of hits from cache equals number from disk.
    - Cache retrieval time is less than disk retrieval time (performance assertion).
  - Notes:
    - This test is marked `#[ignore]` and requires the `LNC_SEEKER_TEST_BAM` environment variable to point at a real BAM file. Run with:

      ```bash
      cd lnc_seeker_lib
      cargo test test_real_bam_cache -- --ignored --nocapture
      ```

## How to run the test-suite

- Unit tests (fast):

  ```bash
  cd lnc_seeker_lib
  cargo test
  ```

- Run a single test by name (use `--test-threads=1` to avoid interference in shared-cache tests):

  ```bash
  cargo test test_cache_integrity -- --test-threads=1
  ```

- Run the ignored integration test (requires `LNC_SEEKER_TEST_BAM` to point at a readable BAM file):

  Set `LNC_SEEKER_TEST_BAM` to a readable BAM path first.

  ```bash
  cargo test test_real_bam_cache -- --ignored --nocapture
  ```

## Dependencies & environment notes

- `test_real_bam_cache` requires `LNC_SEEKER_TEST_BAM` to point at a real BAM file; it uses `noodles::bam` to read records.
- The tests use internal types like `ReadInfo`, `get_read_provider`, and `run_analysis` — changing those APIs can invalidate the expectations recorded here.

If you want, I can also add short doc-comments above each test in the source to keep the intent close to the code.
