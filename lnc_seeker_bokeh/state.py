# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import re
import json
import sys
import threading
import uuid
import time
from . import shared_data
from benchmarking.python.benchmark_manager import get_benchmark_manager

def log_safe(state, msg):
    """Thread-safe logging to the UI. Appends to buffer; UI handles display."""
    print(f"[UI LOG] {msg}")
    with state["lock"]:
        state["log_buffer"].append(msg)
        # Trigger update in the current session only
        for doc, session_cb in list(state["active_docs"].items()):
            try:
                doc.add_next_tick_callback(session_cb)
            except Exception:
                pass

def verify_environment(state, lnc_seeker):
    """Verify lnc_seeker was loaded from the virtualenv (best-effort check)."""
    try:
        venv_dir = os.environ.get("VIRTUAL_ENV") or os.path.dirname(os.path.dirname(sys.executable))
        bb_file = getattr(lnc_seeker, "__file__", None)
        if bb_file:
            if not os.path.abspath(bb_file).startswith(os.path.abspath(venv_dir)):
                msg = f"[ENV CHECK] lnc_seeker loaded from {bb_file} (expected under venv {venv_dir})"
                print(msg)
                with state["lock"]:
                    state["log_buffer"].append(msg)
    except Exception as _e:
        print(f"[ENV CHECK] Failed to verify lnc_seeker path: {_e}")

def load_base_data(state):
    t_start = time.time()
    try:
        if not os.path.exists("config.json"):
            return False
        
        # Load config into a local var first to avoid partial state if file is corrupt
        with open("config.json", 'r') as f:
            cfg = json.load(f)
        
        state["config"] = cfg
        state["config_warnings"] = []
        state["optimization_warnings"] = []
        
        # Initialize benchmarking if enabled
        benchmark = get_benchmark_manager(state)

        # Capture all available BAM paths and map them by their display name
        with state["lock"]:
            bam_config = {}

            # Handle explicit bam_collection_dir if provided
            collection_dir = state["config"]["data_selection"].get("bam_collection_dir")
            if collection_dir:
                collection_dir = os.path.abspath(collection_dir)
                if os.path.isdir(collection_dir):
                    manifest_found = False
                    # Check both names, prioritizing 'directory.json'
                    for manifest_name in ["directory.json", "dictionary.json"]:
                        manifest_path = os.path.join(collection_dir, manifest_name)
                        if os.path.exists(manifest_path):
                            try:
                                with open(manifest_path, 'r') as f:
                                    manifest_data = json.load(f)
                                
                                # Resolve relative paths inside the manifest
                                for _region, entry in manifest_data.items():
                                    if isinstance(entry, dict) and "cohorts" in entry:
                                        valid_cohorts = {}
                                        for _c_name, c_info in entry["cohorts"].items():
                                            p = c_info["path"] if isinstance(c_info, dict) else c_info
                                            if p and not os.path.isabs(p):
                                                abs_p = os.path.join(collection_dir, p)
                                            else:
                                                abs_p = p
                                                
                                            if abs_p and os.path.exists(abs_p):
                                                if isinstance(c_info, dict):
                                                    c_info["path"] = abs_p
                                                    valid_cohorts[_c_name] = c_info
                                                else:
                                                    valid_cohorts[_c_name] = abs_p
                                            else:
                                                msg = f"Skipping missing BAM in manifest: {abs_p} (Cohort: {_c_name})"
                                                log_safe(state, msg)
                                                state["config_warnings"].append(msg)
                                        
                                        entry["cohorts"] = valid_cohorts

                                bam_config.update(manifest_data)
                                log_safe(state, f"Loaded BAM collection manifest from {manifest_path}")
                                manifest_found = True
                                break
                            except Exception as e:
                                log_safe(state, f"Failed to load manifest {manifest_path}: {e}")
                    if not manifest_found:
                        log_safe(state, f"No directory.json or dictionary.json found in {collection_dir}")
                else:
                    log_safe(state, f"BAM collection directory not found: {collection_dir}")
            
            # Record timing for loading manifest and resolving paths
            duration = time.time() - t_start
            benchmark.log_event("Initialization", "load_base_data", duration, 
                               details="Manifest loading and gene-hierarchy resolution")

            state["bam_hierarchy"] = bam_config
            state["all_bam_paths"] = []
            
            clean_hierarchy = {}
            for gene, entry in list(state["bam_hierarchy"].items()):
                # Strictly expect: {"cohorts": {...}, "metadata": {...}}
                if not isinstance(entry, dict) or "cohorts" not in entry:
                    log_safe(state, f"Skipping invalid BAM structure for {gene}: Missing 'cohorts'")
                    continue

                if "metadata" not in entry:
                    entry["metadata"] = {}

                # Filter valid paths and populate global name-to-path mapping
                valid_cohorts = {}
                for cohort_name, val in entry.get("cohorts", {}).items():
                    # Support both legacy string path and metadata object
                    p = val["path"] if isinstance(val, dict) else val
                    if p and os.path.exists(p):
                        if p not in state["all_bam_paths"]:
                            state["all_bam_paths"].append(p)
                        valid_cohorts[cohort_name] = val
                    else:
                        msg = f"Skipping missing BAM: {p} (Cohort: {cohort_name})"
                        log_safe(state, msg)
                        state["config_warnings"].append(msg)
                
                if valid_cohorts:
                    entry["cohorts"] = valid_cohorts
                    clean_hierarchy[gene] = entry
            
            state["bam_hierarchy"] = clean_hierarchy
                        
            # Initially, no samples are loaded to satisfy the requirement
            state["config"]["data_selection"]["bam_paths"] = []

        # Ensure defaults in hierarchical structure
        _ensure_config_defaults(state)
        
        # Validate GTF paths and Assembly Report path
        data_sel = state["config"].get("data_selection", {})
        gtf_paths = data_sel.get("gtf_paths", [])
        valid_gtfs = []
        for g in gtf_paths:
            if os.path.exists(g):
                valid_gtfs.append(g)
                # Performance Check: Look for optimized .lba version
                lba_path = f"{g}.lba"
                if not os.path.exists(lba_path):
                    state["optimization_warnings"].append({
                        "msg": f"GTF '{os.path.basename(g)}' is not optimized. Loading will be slow.",
                        "path": g
                    })
            else:
                msg = f"WARNING: GTF file not found: {g}"
                log_safe(state, msg)
                state["config_warnings"].append(msg)
        
        if valid_gtfs != gtf_paths:
            state["config"]["data_selection"]["gtf_paths"] = valid_gtfs
            
        ass_path = data_sel.get("assembly_report_path")
        if ass_path and not os.path.exists(ass_path):
            msg = f"WARNING: Assembly report not found: {ass_path}"
            log_safe(state, msg)
            state["config_warnings"].append(msg)
            # We don't necessarily remove it, but logging it is good.
            
        return True
    except Exception as e:
        log_safe(state, f"Config load error: {e}")
        return False

def _ensure_config_defaults(state):
    if "coverage_and_junctions_profile" not in state["config"]:
        state["config"]["coverage_and_junctions_profile"] = {}
    if "plot_height" not in state["config"]["coverage_and_junctions_profile"]:
        state["config"]["coverage_and_junctions_profile"]["plot_height"] = 450
    if "min_mapping_quality" not in state["config"]["coverage_and_junctions_profile"]:
        state["config"]["coverage_and_junctions_profile"]["min_mapping_quality"] = 10
    if "adapt_global_y_range_to_normalization" not in state["config"]["coverage_and_junctions_profile"]:
        state["config"]["coverage_and_junctions_profile"]["adapt_global_y_range_to_normalization"] = True
    if "normalize_junctions_by_samples" not in state["config"]["coverage_and_junctions_profile"]:
        state["config"]["coverage_and_junctions_profile"]["normalize_junctions_by_samples"] = True
    
    if "data_selection" not in state["config"]:
        state["config"]["data_selection"] = {}
    if "filter_outliers" not in state["config"]["data_selection"]:
        state["config"]["data_selection"]["filter_outliers"] = False
    if "filter_annotations" not in state["config"]["data_selection"]:
        state["config"]["data_selection"]["filter_annotations"] = False
    
    if "genome_annotations" not in state["config"]:
        state["config"]["genome_annotations"] = {}
    if "plot_height" not in state["config"]["genome_annotations"]:
        state["config"]["genome_annotations"]["plot_height"] = 350
    if "show_intron_direction" not in state["config"]["genome_annotations"]:
        state["config"]["genome_annotations"]["show_intron_direction"] = False
    if "show_full_range" not in state["config"]["genome_annotations"]:
        state["config"]["genome_annotations"]["show_full_range"] = False

    if "full_read_layout" not in state["config"]:
        state["config"]["full_read_layout"] = {}
    if "plot_height" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["plot_height"] = 300

    if "show_reads_line_width" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["show_reads_line_width"] = 2

    if "vertical_gap" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["vertical_gap"] = 4
    
    if "reads_vertical_squeeze" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["reads_vertical_squeeze"] = 0.2
    if "min_vertical_squeeze" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["min_vertical_squeeze"] = 0.05
    if "max_vertical_squeeze" not in state["config"]["full_read_layout"]:
        state["config"]["full_read_layout"]["max_vertical_squeeze"] = 1.2
    
    if "transcript_creator" not in state["config"]:
        state["config"]["transcript_creator"] = {}
    if "show_intron_direction" not in state["config"]["transcript_creator"]:
        state["config"]["transcript_creator"]["show_intron_direction"] = False
    
    return True

def strip_id(name):
    """Removes 'N. ' prefix from sample names for data matching."""
    if not isinstance(name, str):
        return name
    return re.sub(r'^\d+\.\s+', '', name)
