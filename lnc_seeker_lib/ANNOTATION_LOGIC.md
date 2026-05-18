# Annotation Retrieval Logic

This document describes the two-pass retrieval logic implemented in the Rust backend to ensure biological data integrity and prevent "cutted" transcripts in the Genome Annotations diagram.

## Overview

The `get_annotations` function in `lib.rs` uses a two-pass approach when querying indexed GTF files (Tabix/CSI). This ensures that if any part of a transcript overlaps the primary Region of Interest (ROI), the *entire* transcript (including exons far outside the ROI) is retrieved and displayed.

## Two-Pass Implementation

### 1. Pass 1: Discovery (ROI)
The backend first queries the GTF index using the **Core Range**. 
- **Core Range**: Either the metadata-defined "Extracted Range" (default) or the "Full Range" (when checked).
- **Goal**: Identify the `transcript_id` of every record that intersects this range.
- **Output**: A set of unique `transcript_id`s present in the user's primary view.

### 2. Pass 2: Collection (Padded Range)
The backend performs a second query using the **Fetch Range**.
- **Fetch Range**: The Core Range plus a **500kb padding** on both sides.
- **Saturating Logic**: Uses `saturating_sub` to prevent coordinate underflow at the start of chromosomes.
- **Filtering**: Only records belonging to the IDs discovered in Pass 1 are collected.
- **Goal**: Retrieve distal exons and UTRs that belong to the "relevant" transcripts but lie outside the primary ROI.

## Coordinate Handling

- **ROI Focus**: When "Full Range Annotations" is unselected, the system explicitly uses coordinates from `config.json` or the **BAM-derived cache metadata** as the Core Range.
- **GTF Offsets**: Supports per-file coordinate offsets defined in the configuration.
- **Coordinate Consistency**: The system normalizes GTF (1-based) coordinates to match the 0-based system used in BAM files and the lnc-seeker cache.
- **Limit**: A global limit of **25,000 features** is enforced across all processed GTF files to maintain UI performance.

## Fallback
For non-indexed GTF files, the system falls back to a **Streaming Parser**. Due to performance constraints on large files, the streaming parser uses a single-pass filter against the Core Range. For best results (uncut transcripts), using indexed (`.tbi` or `.csi`) GTF files is strongly recommended.
