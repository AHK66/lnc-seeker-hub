# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import gc
import json
import sys
import time
import lnc_seeker
from lnc_seeker_bokeh.state import log_safe
from lnc_seeker_bokeh.data_utils import process_analysis_data
from benchmarking.python.benchmark_manager import get_benchmark_manager


def trim_malloc_heap_linux_only():
    if sys.platform != "linux":
        return

    try:
        import ctypes

        libc = ctypes.CDLL("libc.so.6")
        malloc_trim = libc.malloc_trim
        malloc_trim.argtypes = [ctypes.c_size_t]
        malloc_trim.restype = ctypes.c_int
        malloc_trim(0)
    except Exception:
        pass

def run_analysis_thread(state):
    """
    Main background thread for performing the analysis.
    Uses the Mailbox approach: updates state and notifies UI via callbacks.
    """
    try:
        # Clear any previous fatal errors when a new analysis starts
        with state["lock"]:
            state["fatal_error"] = None

        # Keep track of last known memory values to prevent "zero valleys" in charts
        last_cache_core_mb = 0.0
        last_cache_related_mb = 0.0
        last_cache_annotation_mb = 0.0
        last_json_overhead_mb = 0.0
        last_gene = None

        while True:
            # Check if config is available
            with state["lock"]:
                if state["config"] is None:
                    print("[PIPELINE] Config not loaded yet. Aborting.")
                    state["analysis_running"] = False
                    return
                
                state["analysis_running"] = True
                current_config_json = json.dumps(state["config"])
                mq_val = state["config"]["coverage_and_junctions_profile"].get("min_mapping_quality", 0)
                selected_gtfs = state["config"]["data_selection"].get("selected_gtfs", [])

            log_safe(state, f"Starting Rust analysis [MQ={mq_val}]...")
            benchmark = get_benchmark_manager(state)
            
            # 1. Run Rust Backend
            t_json_start = time.time()
            current_config_json = json.dumps(state["config"])
            t_json_end = time.time()
            serialization_time = t_json_end - t_json_start
            
            # Calculate metrics for benchmarking
            processed_reads = 0
            gene_name = state["config"]["general"].get("gene_name", "")
            
            # If the gene has changed, reset the "sticky" memory values
            if gene_name != last_gene:
                last_cache_core_mb = 0.0
                last_cache_related_mb = 0.0
                last_cache_annotation_mb = 0.0
                last_gene = gene_name

            # Query the core cache size directly from Rust if possible
            try:
                import lnc_seeker
                cache_core_mb = lnc_seeker.get_cache_core_size_py()
                cache_annotation_mb = lnc_seeker.get_cache_annotation_size_py()
            except AttributeError:
                cache_core_mb = last_cache_core_mb
                cache_annotation_mb = last_cache_annotation_mb
            
            # Carry forward the last known related cache (transient data)
            cache_related_mb = last_cache_related_mb

            details = f"Gene: {gene_name}" if gene_name else "System-wide/Refresh Analysis"
            
            selected_samples = state["config"]["data_selection"].get("selected_samples", [])
            hierarchy = state.get("bam_hierarchy", {})
            if gene_name in hierarchy:
                cohorts = hierarchy[gene_name].get("cohorts", {})
                for s in selected_samples:
                    if s in cohorts:
                        c_info = cohorts[s]
                        if isinstance(c_info, dict):
                            processed_reads += c_info.get("num_reads", 0)
                            # Only add if we don't have a direct Rust measurement yet
                            if not cache_core_mb:
                                cache_core_mb += c_info.get("cache_size_bytes", 0) / (1024 * 1024)

            print(f"[BENCHMARK] Python: Config serialization took {serialization_time:.4f}s")
            benchmark.log_event("Processing", "config_serialization", serialization_time, 
                               reads=processed_reads, cache_core_mb=cache_core_mb, 
                               cache_related_mb=cache_related_mb, 
                               cache_annotation_mb=cache_annotation_mb, 
                               json_overhead_mb=last_json_overhead_mb, details=details)

            t0 = time.time()
            # The progress tracker is thread-safe (Arc<ProgressData> inside)
            res_json = lnc_seeker.run_analysis_py(current_config_json, state["progress_tracker"])
            t1 = time.time()
            rust_time = t1 - t0

            # Calculate JSON Overhead in MB (Raw string size)
            json_overhead_mb = len(res_json) / (1024 * 1024)
            last_json_overhead_mb = json_overhead_mb
            with state["lock"]:
                state["last_json_overhead_mb"] = json_overhead_mb

            # Staleness Guard: Check if config changed while Rust was busy
            with state["lock"]:
                if json.dumps(state["config"]) != current_config_json:
                    log_safe(state, "Config changed during Rust analysis. Discarding stale results...")
                    continue

            log_safe(state, f"Rust analysis took {rust_time:.2f}s. Parsing JSON...")
            print(f"[BENCHMARK] Python: Rust backend total (including serialization) took {rust_time:.4f}s")
            
            # 2. Parse Results
            t_parse_start = time.time()
            new_data = json.loads(res_json)
            t_parse_end = time.time()
            deserialization_time = t_parse_end - t_parse_start
            
            # Synchronize with actual cache usage reported from Rust core
            # We do this BEFORE logging 'rust_backend_core' so the report reflects the peak after analysis
            cache_core_mb = new_data.get('cache_core_mb', cache_core_mb)
            cache_related_mb = new_data.get('cache_related_mb', 0.0)
            cache_annotation_mb = new_data.get('cache_annotation_mb', 0.0)
            
            # Update last known values for the next preamble
            last_cache_core_mb = cache_core_mb
            last_cache_related_mb = cache_related_mb
            last_cache_annotation_mb = cache_annotation_mb

            benchmark.log_event("Processing", "rust_backend_core", rust_time, 
                               reads=processed_reads, cache_core_mb=cache_core_mb, 
                               cache_related_mb=cache_related_mb, 
                               cache_annotation_mb=cache_annotation_mb, 
                               json_overhead_mb=json_overhead_mb, details=details)

            print(f"[BENCHMARK] Python: Results deserialization took {deserialization_time:.4f}s")
            benchmark.log_event("Processing", "results_deserialization", deserialization_time, 
                               reads=processed_reads, cache_core_mb=cache_core_mb, 
                               cache_related_mb=cache_related_mb, 
                               cache_annotation_mb=cache_annotation_mb, 
                               json_overhead_mb=json_overhead_mb, details=details)
            
            # 3. Save to state and process
            with state["lock"]:
                state["analysis_data"] = new_data
                state["data_gene_name"] = state["config"]["general"].get("gene_name")
                
            t_proc_start = time.time()
            process_analysis_data(state)
            t_proc_end = time.time()
            processing_time = t_proc_end - t_proc_start
            print(f"[BENCHMARK] Python: process_analysis_data took {processing_time:.4f}s")
            
            # --- MEMORY OPTIMIZATION ---
            # Prune the raw dictionary to save hundreds of MBs of redundant Python objects.
            # We keep 'samples' as it's used by ui_manager for plot updates.
            with state["lock"]:
                if state["analysis_data"]:
                    state["analysis_data"].pop("annotations", None)
            
            # Explicitly clear potentially massive intermediate objects and trigger GC
            del res_json
            del new_data
            gc.collect()
            trim_malloc_heap_linux_only()
            # ---------------------------

            benchmark.log_event("Processing", "process_analysis_data", processing_time, 
                               reads=processed_reads, cache_core_mb=cache_core_mb, 
                               cache_related_mb=cache_related_mb, 
                               cache_annotation_mb=cache_annotation_mb, 
                               json_overhead_mb=json_overhead_mb, details=details)

            # Check for Change Coalescing and set finished state BEFORE notifying
            # This prevents the race condition where UI sees 'running=True' during notification
            with state["lock"]:
                if json.dumps(state["config"]) != current_config_json:
                    log_safe(state, "Configuration changed during analysis. Re-starting...")
                    continue
                else:
                    state["analysis_running"] = False
                    is_final_pass = True

            log_safe(state, "Analysis complete. Notifying UI...")
            
            # 4. Notify UI (The Mailbox delivery)
            with state["lock"]:
                callbacks = list(state["active_docs"].values())
                
            for cb in callbacks:
                try:
                    # Every callback here is bound to a specific document
                    # so this is perfectly isolated.
                    # The closure refresh_session_ui in ui_manager sees its own self.state.
                    cb() 
                except Exception as e:
                    print(f"[PIPELINE] Callback error: {e}")
            
            if is_final_pass:
                break
                    
    except Exception as e:
        import traceback
        err_msg = f"Pipeline Error: {e}\n{traceback.format_exc()}"
        print(err_msg)
        
        # Specific handling for missing files (serious configuration error)
        # On Windows, os error 3 is "The system cannot find the path specified" (Path Not Found)
        err_str = str(e).lower()
        is_path_error = (isinstance(e, OSError) and (getattr(e, 'errno', None) in [2, 3])) or \
                        ("os error 2" in err_str) or ("os error 3" in err_str)
        
        with state["lock"]:
            state["analysis_running"] = False
            # Check for Change Coalescing
            try:
                rerun = json.dumps(state["config"]) != current_config_json
            except NameError:
                rerun = False

            if is_path_error:
                fatal_msg = f"<b>Fatal Configuration Error:</b><br/>{e}<br/><br/>The system cannot find a required file or path. Please verify all BAM and GTF paths in <code>config.json</code> or the linked manifest."
                state["fatal_error"] = fatal_msg
                log_safe(state, "FATAL ERROR: Required file/path missing. Check dashboard.")
            else:
                log_safe(state, f"An error occurred during analysis: {e}")
            
            if rerun:
                log_safe(state, "Config changed during scan. Starting follow-up scan...")
                # The while loop will continue if we didn't break, but we are in except.
                # Since we want to rerun, we should ideally NOT have left the loop, 
                # but for simplicity in this thread-based architecture, we'll let the next trigger start it
                # or we could recursively call, but that's risky. 
                # Actually, the background thread is started by the UI manager.
    finally:
        with state["lock"]:
            state["analysis_running"] = False
