# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import time
import datetime
import threading
import psutil
import functools
from bokeh.plotting import curdoc

class BenchmarkManager:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(BenchmarkManager, cls).__new__(cls)
                cls._instance._initialized = False
            return cls._instance

    def __init__(self, enabled=False, session_id="default"):
        if self._initialized:
            return
        self.enabled = enabled
        self.session_id = session_id
        # Define log file path relative to current repo root
        self.log_file = os.path.join("benchmarking", "data", "performance_log.csv")
        self._write_header()
        self._initialized = True

    def _write_header(self):
        if not self.enabled:
            return
        # If file exists, check if it has the current header
        header = "Timestamp,SessionID,Category,Event,Duration_s,Memory_MB,CacheRelated_MB,CacheCore_MB,CacheAnnotation_MB,JsonOverhead_MB,ProcessedReads,Details\n"
        if not os.path.exists(self.log_file):
            with open(self.log_file, "w") as f:
                f.write(header)
        else:
            # Check for legacy header and update if necessary
            try:
                with open(self.log_file, "r") as f:
                    first_line = f.readline()
                if first_line and "JsonOverhead_MB" not in first_line:
                    # Rename old log and start fresh with new header
                    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                    os.rename(self.log_file, f"performance_log_{timestamp}.csv.bak")
                    with open(self.log_file, "w") as f:
                        f.write(header)
            except Exception:
                pass

    def get_memory_mb(self):
        """Returns current process Resident Set Size in MiB."""
        try:
            process = psutil.Process(os.getpid())
            return process.memory_info().rss / (1024 * 1024)
        except Exception:
            return 0.0

    def log_event(self, category, event, duration=0.0, reads=0, cache_related_mb=0.0, cache_core_mb=0.0, cache_annotation_mb=0.0, json_overhead_mb=0.0, details=""):
        if not self.enabled:
            return
        
        # 1. Fetch Ground Truth from Rust for Caches where possible
        # This prevents "zero valleys" during transitions or interaction events
        if cache_core_mb == 0.0:
            try:
                import lnc_seeker
                cache_core_mb = lnc_seeker.get_cache_core_size_py()
            except (ImportError, AttributeError):
                cache_core_mb = getattr(self, "_last_cache_core_mb", 0.0)
        
        if cache_core_mb > 0:
            self._last_cache_core_mb = cache_core_mb

        if cache_annotation_mb == 0.0:
            try:
                import lnc_seeker
                cache_annotation_mb = lnc_seeker.get_cache_annotation_size_py()
            except (ImportError, AttributeError):
                cache_annotation_mb = getattr(self, "_last_cache_annotation_mb", 0.0)
        
        if cache_annotation_mb > 0:
            self._last_cache_annotation_mb = cache_annotation_mb

        # 2. Related Cache Memory (Coverage/Transient) - Sticky behavior
        if cache_related_mb == 0.0:
            cache_related_mb = getattr(self, "_last_cache_related_mb", 0.0)
        
        if cache_related_mb > 0:
            self._last_cache_related_mb = cache_related_mb

        # 3. JSON Serialization Overhead - Sticky behavior
        if json_overhead_mb == 0.0:
            json_overhead_mb = getattr(self, "_last_json_overhead_mb", 0.0)
        
        if json_overhead_mb > 0:
            self._last_json_overhead_mb = json_overhead_mb

        # Try to resolve session ID from Bokeh context if still "default"
        session = self.session_id
        if session == "default":
            try:
                doc = curdoc()
                if doc and doc.session_context:
                    session = doc.session_context.id
            except Exception:
                pass

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
        memory = self.get_memory_mb()
        
        # Refine empty details with context-aware fallbacks
        if not details:
            if category == "Interaction":
                details = f"UI Action: {event}"
            elif category == "Processing":
                details = "Data processing task"
            elif category == "Initialization":
                details = "System initialization"
            else:
                details = f"Logged {event}"

        # Escape details for CSV
        details_escaped = str(details).replace('"', '""')
        
        line = f"{timestamp},{session},{category},{event},{duration:.4f},{memory:.2f},{cache_related_mb:.2f},{cache_core_mb:.2f},{cache_annotation_mb:.2f},{json_overhead_mb:.2f},{reads},\"{details_escaped}\"\n"
        
        with open(self.log_file, "a") as f:
            f.write(line)

    def wrap_callback(self, category, event_name):
        """Decorator/Wrapper for Bokeh callbacks to automatically time and log them."""
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                if not self.enabled:
                    return func(*args, **kwargs)
                
                t0 = time.time()
                try:
                    result = func(*args, **kwargs)
                    return result
                finally:
                    duration = time.time() - t0
                    # Default metrics
                    reads = 0
                    cache_related_mb = 0.0
                    cache_core_mb = 0.0
                    cache_annotation_mb = 0.0
                    json_overhead_mb = 0.0
                    details = ""
                    
                    # Resolve the target object (usually VisualizerApp instance)
                    target = getattr(func, '__self__', None)
                    if not target and args:
                        target = args[0]
                    
                    # Probe the target for data metrics and context AFTER execution
                    if target and hasattr(target, 'L'):
                        if hasattr(target, 'get_processed_reads'):
                            reads = target.get_processed_reads()
                        if hasattr(target, 'get_cache_related_mb'):
                            cache_related_mb = target.get_cache_related_mb()
                        if hasattr(target, 'get_cache_core_mb'):
                            cache_core_mb = target.get_cache_core_mb()
                        if hasattr(target, 'get_cache_annotation_mb'):
                            cache_annotation_mb = target.get_cache_annotation_mb()
                        if hasattr(target, 'get_json_overhead_mb'):
                            json_overhead_mb = target.get_json_overhead_mb()
                        
                        L = getattr(target, 'L', None)
                        if L and 'sel_gene' in L:
                            gene = L['sel_gene'].value
                            if gene:
                                details = f"Gene: {gene}"
                            else:
                                details = f"Action: {event_name} (No gene selected)"
                        else:
                            details = f"App context: {event_name}"
                    else:
                        details = f"Callback: {event_name}"

                    self.log_event(category, event_name, duration, reads=reads, 
                                   cache_related_mb=cache_related_mb, cache_core_mb=cache_core_mb, 
                                   cache_annotation_mb=cache_annotation_mb, 
                                   json_overhead_mb=json_overhead_mb, details=details)
            return wrapper
        return decorator

def get_benchmark_manager(state=None):
    """Utility to get or initialize the singleton manager from state."""
    enabled = False
    session_id = "default"
    if state:
        session_id = state.get("session_id", "default")
        if "config" in state and state["config"]:
            enabled = state["config"].get("general", {}).get("enable_benchmarking", False)
    
    return BenchmarkManager(enabled=enabled, session_id=session_id)
