# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
# shared_data.py
# This module provides the structure for the session-specific state.
import threading
from collections import deque


def _make_progress_tracker():
    """Create a SessionProgress instance if `lnc_seeker` is available.

    This avoids import-time failures when the compiled extension isn't
    installed (e.g., during development without running `maturin`).
    """
    try:
        import lnc_seeker
        return lnc_seeker.SessionProgress()
    except Exception:
        return None


def create_session_state():
    return {
        # Application configuration loaded from config.json
        "config": None,
        # Raw result data returned from the Rust analysis backend
        "analysis_data": None,
        # The gene name that the current analysis_data corresponds to
        "data_gene_name": None,
        # Tracks progress of long-running analysis tasks (may be None)
        "progress_tracker": _make_progress_tracker(),
        # Set of (start, end) tuples for introns from curated annotations
        "known_introns": set(),
        # Set of (start, end) tuples for all introns in annotations
        "all_introns": set(),
        # Genomic coordinates of slice sites from curated annotations
        "known_sites": set(),
        # Genomic coordinates of slice sites from predicted/uncurated annotations
        "predicted_sites": set(),
        # Boolean flag indicating if an analysis is currently executing
        "analysis_running": False,
        # Circular buffer for log messages displayed in the UI
        "log_buffer": deque(maxlen=500),
        # List of configuration warnings found during path validation
        "config_warnings": [],
        # Map of file stems to cohort names for backend data resolution
        "stem_to_cohort": {},
        # Map of cohort names to absolute BAM file paths
        "cohort_to_path": {},
        # Cache key for avoiding redundant annotation processing
        "last_ann_cache_key": None,
        # If a serious configuration error occurs, this stores the error message
        "fatal_error": None,
        # Map of Bokeh Document objects to their UI refresh callbacks
        "active_docs": {}, # doc -> update_callback
        # Reentrant lock for thread-safe access to this state dictionary
        "lock": threading.RLock(),
        # Processed exon/transcript data for visualization (pandas DataFrame)
        "processed_annotations": None,
        # Flat exon model for core logic and junction assessment (pandas DataFrame)
        "flat_annotations": None,
        # Processed coding sequence regions (pandas DataFrame)
        "processed_cds": None,
        # Metadata and layout coordinates for transcripts (pandas DataFrame)
        "processed_transcripts": None,
        # Labels and positions for gene names in the plot (pandas DataFrame)
        "processed_gene_labels": None,
        # Markers for introns (arrows/directional)
        "processed_markers": None,
        # Intron segments for hover information
        "processed_introns": None,
    }
