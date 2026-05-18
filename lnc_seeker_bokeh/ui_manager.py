# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import sys
import json
import time
import math
import threading
import datetime
import traceback
import numpy as np
import pandas as pd
from bokeh.plotting import figure, curdoc
from bokeh.models import (
    ColumnDataSource, HoverTool, Range1d, LinearColorMapper, ColorBar,
    FixedTicker, Div, Button, Select, MultiSelect, Row, Column, Spinner,
    CheckboxGroup, Slider, LogAxis, BasicTickFormatter, LogTicker,
    LabelSet, NumeralTickFormatter, WheelZoomTool, PanTool, TapTool,
    WheelPanTool, RadioButtonGroup, Toggle, TextInput, NumericInput,
    Span, CustomJS, Tabs, TabPanel, FileInput,
    DataTable, TableColumn, CheckboxEditor, StringEditor, SelectEditor
)
from bokeh.events import Reset, MouseWheel, ButtonClick
from bokeh.layouts import column, row
import lnc_seeker
from lnc_seeker_bokeh.state import log_safe, verify_environment, load_base_data
from lnc_seeker_bokeh.constants import (
    PROGRESS_STYLES, RedGrayBlue11, set_progress_in_progress, 
    set_progress_complete, set_progress_redrawing, get_progress_html,
    set_progress_message, clear_progress
)
from lnc_seeker_bokeh.data_utils import process_analysis_data, get_marked_sets, calculate_global_ranges
from lnc_seeker_bokeh.pipeline import run_analysis_thread
from lnc_seeker_bokeh.plotting_base import get_j_color, get_transcript_color, add_crosshair_to_plot, get_status_style
from lnc_seeker_bokeh.coverage_plot import create_sample_plot, update_sample_data, get_or_update_sample_plot
from lnc_seeker_bokeh.transcript_creator import TranscriptCreator
from lnc_seeker_bokeh.reads_manager import ReadsManager
from lnc_seeker_bokeh.genome_manager import GenomeManager
from lnc_seeker_bokeh.rules_manager import RulesManager
from lnc_seeker_bokeh.selection_manager import SelectionManager
from lnc_seeker_bokeh.ui_layout import initialize_local_state, setup_ui
from benchmarking.python.benchmark_manager import get_benchmark_manager

def _format_mio(n):
    """Formats a number as human-friendly string with Mio/k suffixes."""
    try:
        n = float(n)
        if n >= 1_000_000:
            return f"{n / 1_000_000:.2f} Mio"
        if n >= 1_000:
            return f"{n / 1_000:.1f}k"
        return f"{int(n):,}"
    except (ValueError, TypeError):
        return "N/A"

def _format_size(bytes):
    """Formats bytes as human-friendly size string."""
    try:
        bytes = float(bytes)
        if bytes <= 0:
            return "0 B"
        units = ("B", "KB", "MB", "GB", "TB")
        i = int(math.floor(math.log(bytes, 1024)))
        return f"{bytes / math.pow(1024, i):.2f} {units[i]}"
    except (ValueError, TypeError, OverflowError):
        return "N/A"

def _format_coverage(val):
    """Formats coverage as human-friendly string with 2 decimal places."""
    try:
        v = float(val)
        return f"{v:.2f}x"
    except (ValueError, TypeError):
        return "N/A"

class VisualizerApp:
    def __init__(self, doc):
        self.doc = doc
        self.session_id = doc.session_context.id if doc.session_context else "default"
        
        # In the new explicit architecture, each doc owns its state.
        # If multiple tabs open for the same session, they get independent states
        # to prevent clobbering.
        from . import shared_data
        self.state = shared_data.create_session_state()
        self.state["session_id"] = self.session_id
        
        # Load config immediately so it's available for layout initialization
        load_base_data(self.state)

        # Register document
        with self.state["lock"]:
            def notify_callback():
                try:
                    self.doc.add_next_tick_callback(self.refresh_session_ui)
                except Exception:
                    pass
            self.state["active_docs"][self.doc] = notify_callback

        self.L = initialize_local_state(self)
        self.L["doc"] = doc
        self.L["session_id"] = self.session_id

        # Display configuration warnings if any
        if self.state.get("config_warnings") or self.state.get("optimization_warnings"):
            warn_items = []
            if self.state.get("config_warnings"):
                warn_items.append("<b>File integrity issues:</b>")
                warn_items.extend([f"<li>{w}</li>" for w in self.state["config_warnings"]])
            
            if self.state.get("optimization_warnings"):
                if warn_items: warn_items.append("<br>")
                warn_items.append("<b>Performance Recommendations:</b>")
                for opt in self.state["optimization_warnings"]:
                    msg = opt["msg"]
                    path = opt["path"]
                    warn_items.append(f"<li>{msg} (Run: <code>cargo run -p lnc_seeker_lib -- --optimize-gtf \"{path}\"</code>)</li>")

            warn_html = "<ul>" + "".join(warn_items) + "</ul>"
            self.L["div_config_warnings"].text = f"""
                <div style="margin: 10px; padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 4px; color: #856404; font-family: sans-serif;">
                    <h4 style="margin-top: 0; color: #856404;">⚠️ Configuration & Performance Warnings</h4>
                    {warn_html}
                    <p style="margin-bottom: 0; font-size: 0.9em; font-style: italic;">Optimization is highly recommended for large datasets (e.g., GENCODE).</p>
                </div>
            """
            self.L["div_config_warnings"].visible = True
        
        # Initialize sub-components
        self.creator = TranscriptCreator(self)
        self.reads_manager = ReadsManager(self)
        self.genome_manager = GenomeManager(self)
        self.rules_manager = RulesManager(self)
        self.selection_manager = SelectionManager(self)
        
        # Initialize benchmarking manager for this session
        self.benchmark = get_benchmark_manager(self.state)
        self._last_verified_cache_mb = 0.0
        self._last_verified_related_mb = 0.0

        # Consistent mapping for tools to access managers from L
        self.L["reads_manager"] = self.reads_manager
        self.L["selection_manager"] = self.selection_manager
        
        setup_ui(self)
        self._setup_callbacks()
        self.update_normalization_ui()

        # Handle session cleanup to prevent reference errors to dead documents
        self.doc.on_session_destroyed(self.on_session_destroyed)

    def on_session_destroyed(self, session_context):
        """Cleanup session registrations when the user closes the tab."""
        try:
            with self.state["lock"]:
                if self.doc in self.state["active_docs"]:
                    del self.state["active_docs"][self.doc]
            print(f"[SESSION] Cleaned up state for session: {self.session_id}")
        except Exception:
            pass
        
    def on_reset_zoom_click(self):
        """Resets x-range to the core area (min_x..max_x) detected by Rust and y-range to current max depth."""
        with self.state["lock"]:
            data = self.state.get("analysis_data", {})
            min_x = data.get("min_x")
            max_x = data.get("max_x")
            
            if min_x is not None and max_x is not None and max_x > min_x:
                log_safe(self.state, f"Resetting view to core area: {min_x}..{max_x}")
                self.L["shared_x_range"].start = min_x
                self.L["shared_x_range"].end = max_x

            # Reset Y ranges for all coverage plots to the current locked_y_top
            y_max_list = self.L["ds_core_range"].data.get('y_max', [])
            if y_max_list:
                y_target = y_max_list[0]
                for name, res in self.L["sample_plots"].items():
                    p = res[0]
                    p.y_range.start = 0.1
                    p.y_range.end = y_target

    def on_reset_ann_click(self):
        self.genome_manager.on_reset_ann_click()

    def schedule_update_all_samples(self, set_rendered=False, log_msg=None):
        """Coalesce and run plot refreshes on the document thread."""
        if self.L.get("update_pending", False):
            return
        self.L["update_pending"] = True

        def _run():
            # Guard against updates on a document that is being destroyed or already gone
            if not self.doc or self.doc.session_context is None:
                self.L["update_pending"] = False
                return
                
            try:
                self.update_all_samples()
                if set_rendered:
                    self.L["data_rendered"] = True
                    if log_msg:
                        log_safe(self.state, log_msg)
            except Exception as e:
                import traceback
                print(f"Error in update_all_samples: {e}\n{traceback.format_exc()}")
            finally:
                self.L["update_pending"] = False

        self.doc.add_next_tick_callback(_run)

    def refresh_session_ui(self):
        """Session-safe UI refresh for logs and status."""
        if not self.doc or self.doc.session_context is None:
            return
            
        try:
            # Check for fatal errors first
            with self.state["lock"]:
                fatal_err = self.state.get("fatal_error")
                finished = (not self.state.get("analysis_running", False) and self.state["analysis_data"] is not None)
            
            if fatal_err:
                self.L["div_fatal_error"].text = f"""
                    <div style="padding: 30px; border: 2px solid #dc3545; border-radius: 8px; background-color: #fff8f8; color: #721c24; margin: 20px auto; max-width: 800px; font-family: sans-serif;">
                        <h2 style="color: #c82333; margin-top: 0; display: flex; align-items: center;">
                            <span style="font-size: 1.5em; margin-right: 10px;">⚠️</span> Fatal Configuration Error
                        </h2>
                        <div style="font-size: 1.1em; line-height: 1.5;">
                            {fatal_err}
                        </div>
                        <div style="margin-top: 25px; padding-top: 15px; border-top: 1px solid #f5c6cb; font-size: 0.9em; color: #842029;">
                            <b>Next Steps:</b> Update the configuration file and trigger a re-analysis by selecting a gene or changing a setting.
                        </div>
                    </div>
                """
                self.L["div_fatal_error"].visible = True
                self.L["plot_column"].children = [self.L["div_manual"], self.L["div_config_warnings"], self.L["div_fatal_error"]]
                return

            if finished and not self.L.get("data_rendered", False):
                self.schedule_update_all_samples(set_rendered=True, log_msg="Plots updated with final analysis results.")
        except Exception as e:
            print(f"Error in refresh_session_ui: {e}")

    def update_progress(self):
        """Periodic callback (500ms) to update progress bar and detect completion."""
        # sys.stdout.write(".") ; sys.stdout.flush() 
        try:
            # Check if analysis just finished to trigger automatic plot update
            with self.state["lock"]:
                finished = (not self.state.get("analysis_running", False) and self.state["analysis_data"] is not None)
            
            if finished and not self.L.get("data_rendered", False):
                self.schedule_update_all_samples(set_rendered=True, log_msg="Plots refreshed automatically.")

            # Only show progress if THIS session is actually running analysis
            if not self.state.get("analysis_running", False):
                if (self.L.get("is_redrawing", False) or 
                    self.L.get("is_squeezing", False) or 
                    self.L.get("is_fetching_reads", False) or
                    self.L.get("is_sticky_message", False)):
                     # Keep current notification if a redraw, squeeze, fetch or sticky msg is active
                     return
                if self.state["analysis_data"] is not None:
                     set_progress_complete(self.L["div_progress"])
                else:
                     clear_progress(self.L["div_progress"])
                return

            stage, current, total = self.state["progress_tracker"].get_status()
            stage_map = {0: "Idle", 1: "Analyzing BAMs", 2: "Processing GTF", 3: "Finalizing JSON"}
            stage_text = stage_map.get(stage, "Processing")
            
            if stage == 0:
                 set_progress_in_progress(self.L["div_progress"])
                 return

            percent = 0
            if total > 0:
                percent = min(100, int((current / total) * 100))
            
            # Simple HTML progress bar
            bar_color = "#e67e22" if stage == 1 else "#3498db"
            prog_str = f"{percent}%" if total > 0 else f"{current:,} rec"
            text_color = "white" if (percent > 50 or total == 0) else "#495057"
            
            self.L["div_progress"].text = get_progress_html(stage_text, percent if total > 0 else 100, prog_str, bar_color, text_color)
        except Exception as e:
            pass

    def _get_sorted_cohort_data(self, gene_name):
        """Returns a sorted list of cohort metadata for consistent UI display."""
        if not gene_name or gene_name not in self.state.get("bam_hierarchy", {}):
            return []
            
        entry = self.state["bam_hierarchy"][gene_name]
        cohorts = entry.get("cohorts", {})
        
        cohort_data = []
        for name, info in cohorts.items():
            status = info.get("status", "N/A") if isinstance(info, dict) else "N/A"
            tissue = info.get("tissue", "N/A") if isinstance(info, dict) else "N/A"
            num_samples = info.get("num_samples", 1) if isinstance(info, dict) else 1
            num_reads = info.get("num_reads", 0) if isinstance(info, dict) else 0
            avg_coverage = info.get("avg_coverage_per_sample", 0.0) if isinstance(info, dict) else 0.0
            path = info.get("path", info) if isinstance(info, dict) else info
            cohort_data.append({
                "name": name, 
                "status": status, 
                "tissue": tissue, 
                "num_samples": num_samples, 
                "num_reads": num_reads,
                "avg_coverage": avg_coverage,
                "path": path
            })

        # Sort: Status (DESC), then Tissue (ASC), then Name (ASC)
        cohort_data.sort(key=lambda x: x["name"])
        cohort_data.sort(key=lambda x: x["tissue"])
        cohort_data.sort(key=lambda x: x["status"], reverse=True)
        
        return cohort_data

    def _get_gene_info_html(self, gene_name, preselected_samples=None):
        """Generates a comprehensive HTML overview for a selected gene and its cohorts."""
        if not gene_name or gene_name not in self.state.get("bam_hierarchy", {}):
            return ""
        
        preselected_samples = preselected_samples or []
        
        entry = self.state["bam_hierarchy"][gene_name]
        meta = entry.get("metadata", {})
        cohort_data = self._get_sorted_cohort_data(gene_name)
        
        region = meta.get("region", "Unknown Region")
        offset = meta.get("offset", 0)
        
        html = f"""
        <div style="padding: 20px; font-family: sans-serif; color: #2c3e50;">
            <div style="border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-bottom: 20px;">
                <h1 style="margin: 0; color: #2c3e50;">Gene: {gene_name}</h1>
                <p style="margin: 5px 0 0 0; color: #7f8c8d; font-size: 1.1em;">
                    <b>Genomic context:</b> {region} (Offset: {offset} bp)
                </p>
            </div>
            
            <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 20px;">
                <h3 style="color: #2c3e50; margin: 0;">Available Sample Cohorts</h3>
                <button onclick="window.selectFromPanel(this)" style="background-color: #3498db; color: white; border: none; padding: 10px 20px; border-radius: 4px; cursor: pointer; font-weight: 600; transition: background 0.2s;">
                    Select Checked Samples
                </button>
            </div>
            <p style="color: #7f8c8d; margin-top: 10px;">The following samples are available for analysis. Use the checkboxes to select multiple at once, then click "Select Checked Samples".</p>
            
            <table style="width: 100%; border-collapse: collapse; margin-top: 15px; background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border-radius: 4px; overflow: hidden;">
                <thead>
                    <tr style="background-color: #f8f9fa; border-bottom: 2px solid #dee2e6; text-align: left;">
                        <th style="padding: 12px; font-weight: 600; width: 40px; text-align: center;">#</th>
                        <th style="padding: 12px; font-weight: 600; width: 40px; text-align: center;"><input type="checkbox" onclick="let cs=this.closest('table').querySelectorAll('.cohort-chk'); cs.forEach(c=>c.checked=this.checked)"></th>
                        <th style="padding: 12px; font-weight: 600;">Cohort Name</th>
                        <th style="padding: 12px; font-weight: 600;">Tissue</th>
                        <th style="padding: 12px; font-weight: 600;">Status</th>
                        <th style="padding: 12px; font-weight: 600;">Samples</th>
                        <th style="padding: 12px; font-weight: 600;">Read Count</th>
                        <th style="padding: 12px; font-weight: 600;">Coverage per sample</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for i, c in enumerate(cohort_data):
            name = c["name"]
            tissue = c["tissue"]
            status = c["status"]
            num_samples = c["num_samples"]
            num_reads = c["num_reads"]
            avg_coverage = c["avg_coverage"]
            path = c["path"]
            
            # This matches the value format in the MultiSelect: "1. Name"
            val = f"{i+1}. {name}"
            checked_attr = "checked" if val in preselected_samples else ""
            
            # Format numbers for display
            display_reads = _format_mio(num_reads)
            display_coverage = _format_coverage(avg_coverage)
            
            # Get status colors from config
            sc = get_status_style(self.state, status)
            bg_color = sc.get("bg", "#e1f5fe")
            fg_color = sc.get("fg", "#01579b")
            
            html += f"""
                    <tr style="border-bottom: 1px solid #eee;">
                        <td style="padding: 12px; text-align: center; color: #95a5a6; font-size: 0.9em;">{i+1}</td>
                        <td style="padding: 12px; text-align: center;"><input type="checkbox" class="cohort-chk" data-val="{val}" {checked_attr}></td>
                        <td style="padding: 12px; font-weight: 500;" title="Path: {path}">{name}</td>
                        <td style="padding: 12px;">{tissue}</td>
                        <td style="padding: 12px;"><span style="background: {bg_color}; color: {fg_color}; padding: 2px 8px; border-radius: 12px; font-size: 0.85em;">{status}</span></td>
                        <td style="padding: 12px;">{num_samples}</td>
                        <td style="padding: 12px; font-size: 0.9em; color: #34495e;">{display_reads}</td>
                        <td style="padding: 12px; font-size: 0.9em; color: #34495e;">{display_coverage}</td>
                    </tr>
            """
            
        html += """
                </tbody>
            </table>
            
            <div style="margin-top: 30px; padding: 12px 20px; background-color: #fcf8e3; border-left: 5px solid #f0ad4e; color: #8a6d3b;">
                <p style="margin: 0; font-size: 0.95em;"><b>Tip:</b> You can also use the sidebar on the left for quick adjustments without returning to this panel.</p>
            </div>
        </div>
        """
        return html

    def get_processed_reads(self):
        """Calculates total reads for the current selection for benchmarking."""
        from lnc_seeker_bokeh.state import strip_id
        try:
            gene = self.L["sel_gene"].value
            selected = self.L["sel_samples"].value
            if not gene or not selected:
                return 0
            
            hierarchy = self.state.get("bam_hierarchy", {})
            if gene not in hierarchy:
                return 0
            
            cohorts = hierarchy[gene].get("cohorts", {})
            total_reads = 0
            for s in selected:
                raw_name = strip_id(s)
                if raw_name in cohorts:
                    c_info = cohorts[raw_name]
                    if isinstance(c_info, dict):
                        total_reads += c_info.get("num_reads", 0)
            return total_reads
        except Exception:
            return 0

    def get_cache_core_mb(self):
        """Calculates total LNC1 cache size in MiB for the current selection."""
        from lnc_seeker_bokeh.state import strip_id
        try:
            # 1. First priority: Check analysis result if available
            with self.state["lock"]:
                if self.state.get("analysis_data"):
                    # Only use analysis data if it matches the current UI selection
                    if self.state.get("data_gene_name") == self.L["sel_gene"].value:
                        val = self.state["analysis_data"].get("cache_core_mb", 0.0)
                        return val

            # 2. Second priority: Query Rust directly if we are in a transition
            try:
                import lnc_seeker
                return lnc_seeker.get_cache_core_size_py()
            except (ImportError, AttributeError):
                pass

        except Exception:
            pass
        return 0.0

    def get_cache_annotation_mb(self):
        """Calculates annotation cache size (LBA/GTF) in MiB."""
        try:
            # 1. First priority: Check analysis result
            with self.state["lock"]:
                if self.state.get("analysis_data"):
                    if self.state.get("data_gene_name") == self.L["sel_gene"].value:
                        return self.state["analysis_data"].get("cache_annotation_mb", 0.0)

            # 2. Second priority: Query Rust directly
            try:
                import lnc_seeker
                return lnc_seeker.get_cache_annotation_size_py()
            except (ImportError, AttributeError):
                pass
        except Exception:
            pass
        return 0.0

    def get_json_overhead_mb(self):
        """Returns the raw size of the last JSON payload received from Rust."""
        try:
            with self.state["lock"]:
                return self.state.get("last_json_overhead_mb", 0.0)
        except Exception:
            return 0.0

    def get_cache_related_mb(self):
        """Returns the memory overhead for coverage vectors (reported from Rust)."""
        try:
            with self.state["lock"]:
                if self.state.get("analysis_data"):
                    # Only use analysis data if it matches the current UI selection
                    if self.state.get("data_gene_name") == self.L["sel_gene"].value:
                        val = self.state["analysis_data"].get("cache_related_mb", 0.0)
                        self._last_verified_related_mb = val
                        return val
            # Return last known verified value to prevent "zero valleys" during transitions
            return self._last_verified_related_mb
        except Exception:
            return 0.0

    def update_all_samples(self):
        """Main entry point for refreshing all sample plots based on current view/selections."""
        if self.state.get("config") is None or self.state.get("fatal_error"):
            return
        
        t0 = time.time()
        try:
            # Verify that loaded data matches selected gene to prevent race condition mixups
            with self.state["lock"]:
                data_gene = self.state.get("data_gene_name")
                current_gene = self.L["sel_gene"].value
                if data_gene and current_gene and data_gene != current_gene:
                    # If we have data but it's for the wrong gene, it's stale. 
                    # Don't render experimental plots yet.
                    return

            # If no gene is selected, show the manual/intro and return early.
            if self.L["sel_gene"].value == "":
                self.L["div_gene_info"].visible = False
                self.L["plot_column"].children = [self.L["div_manual"], self.L["div_config_warnings"]]
                self.L["last_ann_state"] = (self.L["sel_gene"].value, None, None, None, None, None)
                self.L["btn_show_cohort_selection"].disabled = True
                self.L["btn_show_cohort_selection"].button_type = "default"
                return

            # Prevent accidental clearing of selection during mass data updates
            self.L["selection_updating"] = True
            selected = self.L["sel_samples"].value
            # Update "Show Cohort Selection" button state
            # Active only if plots are shown (panel hidden) AND a gene is selected
            if not selected or self.L.get("force_show_cohort_selection"):
                self.L["btn_show_cohort_selection"].disabled = True
                self.L["btn_show_cohort_selection"].button_type = "default"
            else:
                self.L["btn_show_cohort_selection"].disabled = False
                self.L["btn_show_cohort_selection"].button_type = "primary"

            # --- Parameter State Tracking to avoid redundant data updates ---
            # Annotation state factors
            ann_state = (
                self.L["sel_gene"].value,
                id(self.state.get("processed_annotations")),
                id(self.state.get("processed_cds")),
                id(self.state.get("processed_transcripts")),
                id(self.state.get("processed_gene_labels")),
                id(self.state.get("processed_markers")),
                tuple(self.L["chk_unsupported_introns"].active) if "chk_unsupported_introns" in self.L else (),
                self.state.get("analysis_running", False)
            )
            ann_changed = (ann_state != self.L.get("last_ann_state"))
            self.L["last_ann_state"] = ann_state

            # Sample data state factors
            # selected = self.selection_manager.get_selected_samples()
            
            rules_data = self.L["src_shared_rules"].data
            mark_reqs = {}
            for r_type in ["curated", "predicted", "novel"]:
                presence = [rules_data['sample'][i] for i, r in enumerate(rules_data.get(r_type, [])) if r == "+ (Present)"]
                absence = [rules_data['sample'][i] for i, r in enumerate(rules_data.get(r_type, [])) if r == "- (Absent)"]
                mark_reqs[r_type] = {"presence": presence, "absence": absence}
            
            min_reads = self.L["sld_min_reads"].value
            filter_flanks = (0 in self.L["chk_filter_flanks"].active) if "chk_filter_flanks" in self.L else False
            active_types_idx = tuple(sorted(self.L["chk_types"].active)) if "chk_types" in self.L else ()
            amb_active = tuple(self.L["chk_ambiguity"].active) if "chk_ambiguity" in self.L else ()
            show_full_cov = (0 in self.L["chk_full_cov"].active) if "chk_full_cov" in self.L else False
            show_bg = (0 in self.L["chk_show_bg"].active) if "chk_show_bg" in self.L else False
            normalize = (0 in self.L["chk_normalize"].active) if "chk_normalize" in self.L else False
            normalize_js = (0 in self.L["chk_normalize_junctions"].active) if "chk_normalize_junctions" in self.L else True
            adapt_global_y = (0 in self.L["chk_global_y_range_normalize"].active) if "chk_global_y_range_normalize" in self.L else True
            amb_factor = self.L["sld_ambiguity_factor"].value if "sld_ambiguity_factor" in self.L else 1.0
            
            data_state = (
                tuple(selected),
                tuple(mark_reqs["curated"]["presence"]),
                tuple(mark_reqs["curated"]["absence"]),
                tuple(mark_reqs["predicted"]["presence"]),
                tuple(mark_reqs["predicted"]["absence"]),
                tuple(mark_reqs["novel"]["presence"]),
                tuple(mark_reqs["novel"]["absence"]),
                min_reads,
                filter_flanks,
                active_types_idx,
                amb_active,
                show_full_cov,
                show_bg,
                normalize,
                normalize_js,
                adapt_global_y,
                amb_factor,
                id(self.state.get("analysis_data")),
                self.state.get("analysis_running", False)
            )
            data_changed = (data_state != self.L.get("last_data_state"))
            self.L["last_data_state"] = data_state
            self.L["last_show_full_cov"] = show_full_cov
            self.L["last_show_bg"] = show_bg
            self.L["last_normalize"] = normalize

            # Creator state factors
            cur_junctions = tuple(self.L.get("user_junctions", []))
            creator_state = (
                cur_junctions,
                self.L["num_t_start"].value if "num_t_start" in self.L else None,
                self.L["num_t_end"].value if "num_t_end" in self.L else None,
                self.L["sel_strand"].value if "sel_strand" in self.L else "+",
                self.L["txt_transcript_id"].value if "txt_transcript_id" in self.L else ""
            )
            creator_changed = (creator_state != self.L.get("last_creator_state"))
            
            fixed_width = (0 in self.L["chk_fixed_width"].active) if self.L.get("chk_fixed_width") else False
            target_width = self.L["spn_plot_width"].value if self.L.get("spn_plot_width") else 1200
            cov_target_h = self.L["spn_height"].value if self.L.get("spn_height") else self.state["config"]["coverage_and_junctions_profile"].get("plot_height", 450)
            ann_target_h = self.L["spn_ann_height"].value if self.L.get("spn_ann_height") else self.state["config"]["genome_annotations"].get("plot_height", 350)
            reads_target_h = self.L["spn_reads_height"].value if self.L.get("spn_reads_height") else self.state["config"]["full_read_layout"].get("plot_height", 300)
            backend = self.state["config"]["general"].get("output_backend", "canvas")

            # Update Annotations
            self.genome_manager.update_genome_data(ann_changed, data_changed)

            # Update Extension Metadata Lines
            if ann_changed:
                ext_x, ext_desc, ext_region, ext_color, ext_dash = [], [], [], [], []
                gene = self.L["sel_gene"].value
                if gene in self.state.get("bam_hierarchy", {}):
                    meta = self.state["bam_hierarchy"][gene].get("metadata", {})
                    offset = meta.get("offset", 0)
                    region_str = meta.get("region", "")
                    
                    eb_cfg = self.state.get("config", {}).get("general", {}).get("extension_boundaries", {})
                    core_cfg = eb_cfg.get("core", {"color": "firebrick", "line_width": 2, "line_dash": "solid", "alpha": 0.8})
                    ext_cfg = eb_cfg.get("extraction", {"color": "navy", "line_width": 2, "line_dash": "dashed", "alpha": 0.8})
                    
                    if offset > 0 and ":" in region_str and "-" in region_str:
                        try:
                            _chrom, coords = region_str.split(":")
                            r_start, r_end = map(int, coords.split("-"))
                            # Inner core boundaries
                            ext_x.extend([r_start + offset, r_end - offset])
                            ext_desc.extend(["Core Start", "Core End"])
                            ext_color.extend([core_cfg.get("color", "firebrick")] * 2)
                            ext_dash.extend([core_cfg.get("line_dash", "solid")] * 2)
                            
                            # Outer extraction boundaries
                            ext_x.extend([r_start, r_end])
                            ext_desc.extend(["Extraction Start", "Extraction End"])
                            ext_color.extend([ext_cfg.get("color", "navy")] * 2)
                            ext_dash.extend([ext_cfg.get("line_dash", "dashed")] * 2)
                            
                            ext_region = [region_str] * 4
                            log_safe(self.state, f"Detected extraction boundaries for {gene}: Core={ext_x[:2]}, Full={ext_x[2:]}")
                        except Exception as e:
                            log_safe(self.state, f"Error parsing extension for {gene}: {e}")
                    else:
                        log_safe(self.state, f"No extension/offset metadata found for {gene} (Offset: {offset}, Region: {region_str})")
                else:
                    log_safe(self.state, f"Gene {gene} not found in bam_hierarchy for extension check.")
                
                # Fetch visual alpha/width from config
                eb_cfg = self.state.get("config", {}).get("general", {}).get("extension_boundaries", {})
                w1 = eb_cfg.get("core", {}).get("line_width", 2)
                w2 = eb_cfg.get("extraction", {}).get("line_width", 2)
                a1 = eb_cfg.get("core", {}).get("alpha", 0.8)
                a2 = eb_cfg.get("extraction", {}).get("alpha", 0.8)
                
                ext_widths = [w1, w1, w2, w2] if len(ext_x) == 4 else []
                ext_alphas = [a1, a1, a2, a2] if len(ext_x) == 4 else []

                # Use safe bounds for both log and linear plots; x0/x1 explicitly named
                self.L["ds_extension"].data = dict(
                    x0=ext_x, x1=ext_x, y0=[0.01]*len(ext_x), y1=[1000000]*len(ext_x), 
                    desc=ext_desc, region=ext_region, color=ext_color, dash=ext_dash,
                    width=ext_widths, alpha=ext_alphas
                )
                if ext_x:
                    log_safe(self.state, f"Updated ds_extension with {len(ext_x)} lines.")
                else:
                    log_safe(self.state, f"ds_extension cleared.")

            _active_idx = self.L["chk_types"].active
            show_marked = (0 in _active_idx)
            active_types = [["curated", "predicted", "novel"][i-1] for i in _active_idx if i > 0]
            
            show_markers = (0 in self.L["chk_markers"].active) if "chk_markers" in self.L else False
            show_cds = (0 in self.L["chk_cds"].active) if "chk_cds" in self.L else True
            show_full_cov = (0 in self.L["chk_full_cov"].active) if "chk_full_cov" in self.L else False
            show_cliffs = (0 in self.L["chk_cliffs"].active) if "chk_cliffs" in self.L else True

            data = self.state["analysis_data"] or {}
            analysis_samples = data.get("samples", {})
            marked_sets = get_marked_sets(self.state, analysis_samples, mark_reqs, min_reads)
            
            if data_changed or ann_changed or creator_changed:
                log_safe(self.state, f"Updating layout for selected samples: {selected}")
            
            new_plots_content = []
            if self.L["div_config_warnings"].visible:
                new_plots_content.append(self.L["div_config_warnings"])

            active_tab = self.L["tabs"].active if self.L.get("tabs") else 0
            
            # SHOW INFO PANEL IF NO SAMPLES SELECTED OR FORCED BY BUTTON
            if not selected or self.L.get("force_show_cohort_selection"):
                self.L["div_gene_info"].visible = True
                new_plots_content.append(self.L["div_gene_info"])
                self.L["plot_column"].children = new_plots_content
                return

            self.L["div_gene_info"].visible = False
            
            first_p = None
            normalize = (0 in self.L["chk_normalize"].active) if "chk_normalize" in self.L else False
            adapt_global_y = (0 in self.L["chk_global_y_range_normalize"].active) if "chk_global_y_range_normalize" in self.L else True
            
            gene = self.L["sel_gene"].value if "sel_gene" in self.L else None
            cohort_metadata = {}
            if gene and gene in self.state.get("bam_hierarchy", {}):
                cohort_metadata = self.state["bam_hierarchy"][gene].get("cohorts", {})

            # Recalculate global ranges, optionally adapting to normalization
            global_x_min, global_x_max, global_y_max = calculate_global_ranges(
                data, selected, 
                normalize=(normalize and adapt_global_y), 
                cohort_metadata=cohort_metadata
            )

            locked_y_top = global_y_max * 1.5
            if global_x_min is not None and global_x_max is not None:
                self.L["ds_core_range"].data = dict(start=[global_x_min], end=[global_x_max], y_max=[locked_y_top])
                if self.L.get("reset_x_range", False) or (self.L["shared_x_range"].start == 0 and self.L["shared_x_range"].end == 1):
                    self.L["shared_x_range"].start = global_x_min
                    self.L["shared_x_range"].end = global_x_max
                    self.L["reset_x_range"] = False
            
            for name in selected:
                s_data = self.selection_manager.get_sample_data(name)
                res, is_new_plot = get_or_update_sample_plot(self.state, 
                    name, self.L, self.on_reset_zoom_click, self.on_junction_selection_change,
                    cov_target_h, backend, fixed_width, target_width, locked_y_top,
                    show_full_cov, show_bg, show_cliffs
                )
                p = res[0]
                if not first_p: first_p = p
                
                if data_changed or is_new_plot:
                    if s_data:
                        update_sample_data(self.state, name, res, s_data, self.L, min_reads, active_types, show_marked, marked_sets)
                new_plots_content.append(p)
            
            # Transcript Creator Preview Plot
            if active_tab == 1:
                p_creator = self.creator.get_or_create_preview_plot(fixed_width, target_width, backend)
                self.creator.update_preview_data(creator_state, creator_changed)
                new_plots_content.append(p_creator)

            # Annotation Plot
            p_ann = self.genome_manager.get_or_create_genome_plot(ann_target_h, fixed_width, target_width, backend)
            new_plots_content.append(p_ann)

            p_reads = self.reads_manager.get_or_create_reads_plot(reads_target_h, fixed_width, target_width, backend)
            new_plots_content.append(p_reads)
            restore_jr = self.L.get("last_selected_junctions")
            if restore_jr:
                for s_n, p_os in self.L["sample_plots"].items():
                    ds_js_r = p_os[2]
                    data_r = ds_js_r.data
                    inds_r = []
                    for i in range(len(data_r.get('x0', []))):
                        if (data_r['reference'][i], data_r['x0'][i], data_r['x1'][i]) in restore_jr:
                            inds_r.append(i)
                    if inds_r:
                        self.L["selection_updating"] = True
                        ds_js_r.selected.indices = inds_r
                        break
            self.L["plot_column"].children = new_plots_content
        finally:
            self.L["is_redrawing"] = False
            self.L["selection_updating"] = False
            
            # Explicitly provide context-rich details for the plot refresh
            gene = self.L.get("sel_gene").value if "sel_gene" in self.L else "N/A"
            samples_count = len(self.L.get("sel_samples").value) if "sel_samples" in self.L else 0
            refresh_details = f"Gene: {gene} | Samples: {samples_count}"
            
            # Fix: Pass metrics explicitly to log_event as it's not handled by the decorator here
            reads = self.get_processed_reads()
            c_related = self.get_cache_related_mb()
            c_core = self.get_cache_core_mb()
            c_annotation = self.get_cache_annotation_mb()
            c_json = self.get_json_overhead_mb()
            self.benchmark.log_event("UI", "update_all_samples", time.time() - t0, 
                                   reads=reads, cache_related_mb=c_related, 
                                   cache_core_mb=c_core, cache_annotation_mb=c_annotation, 
                                   json_overhead_mb=c_json, details=refresh_details)

    def on_junction_selection_change(self, source, attr=None, old=None, new=None):
        self.selection_manager.on_junction_selection_change(source, attr, old, new)

    def on_junction_types_change(self, attr, old, new):
        self.L["is_redrawing"] = True
        set_progress_redrawing(self.L["div_progress"])
        self.doc.add_next_tick_callback(self.update_all_samples)

    def on_backend_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        backend = "svg" if "svg" in new.lower() else "canvas"
        if "general" not in self.state["config"]: self.state["config"]["general"] = {}
        self.state["config"]["general"]["output_backend"] = backend
        log_safe(self.state, f"Rendering backend requested: {backend.upper()}")
        self.L["sample_plots"] = {}; self.L["p_ann"] = None; self.L["p_reads"] = None
        self.update_all_samples()

    def on_fixed_width_toggle(self, attr, old, new):
        self.L["sample_plots"] = {}; self.L["p_ann"] = None; self.L["p_reads"] = None
        self.update_all_samples()

    def on_range_change(self, attr, old, new):
        """Handle zoom/pan events to trigger coverage resampling."""
        if self.L.get("selection_updating") or self.L.get("is_redrawing"):
            return

        # Debounce the resampling to avoid lag during smooth panning
        curr_doc = self.doc
        if self.L.get("debounce_resample"):
            try:
                curr_doc.remove_timeout_callback(self.L["debounce_resample"])
            except Exception:
                pass

        def do_resample():
            self.L["debounce_resample"] = None
            from lnc_seeker_bokeh.coverage_plot import update_sample_coverage

            selected = self.L["sel_samples"].value
            for name in selected:
                if name in self.L["sample_plots"]:
                    res = self.L["sample_plots"][name]
                    s_data = self.selection_manager.get_sample_data(name)
                    if s_data:
                        update_sample_coverage(self.state, res, s_data, self.L, name=name)

        self.L["debounce_resample"] = curr_doc.add_timeout_callback(do_resample, 200)

    def on_mq_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        new_mq = int(new)
        if self.state["config"]["coverage_and_junctions_profile"].get("min_mapping_quality") == new_mq: return
        self.state["config"]["coverage_and_junctions_profile"]["min_mapping_quality"] = new_mq
        # Also use the same MQ for ambiguity highlighting
        self.state["config"]["coverage_and_junctions_profile"]["high_ambiguity_highlighting"]["ambiguity_min_mapping_quality"] = new_mq
        
        curr_doc = self.doc
        if self.L.get("debounce_mq"):
            try: curr_doc.remove_timeout_callback(self.L["debounce_mq"])
            except: pass
        def trigger_analysis():
            self.L["debounce_mq"] = None; should_start_thread = False
            with self.state["lock"]:
                self.state["analysis_data"] = None
                self.state["data_gene_name"] = None
                if self.state.get("analysis_running", False): self.L["data_rendered"] = False
                else: self.L["data_rendered"] = False; self.state["analysis_running"] = True; should_start_thread = True; set_progress_in_progress(self.L["div_progress"])
            if should_start_thread: threading.Thread(target=run_analysis_thread, args=(self.state,), daemon=True).start()
        self.L["debounce_mq"] = curr_doc.add_timeout_callback(trigger_analysis, 500)

    def on_ambiguity_factor_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        new_val = float(new)
        if self.state["config"]["coverage_and_junctions_profile"]["high_ambiguity_highlighting"]["ambiguity_highlight"].get("threshold") == new_val:
            return
        self.state["config"]["coverage_and_junctions_profile"]["high_ambiguity_highlighting"]["ambiguity_highlight"]["threshold"] = new_val
        
        curr_doc = self.doc
        if self.L.get("debounce_amb_factor"):
            try: curr_doc.remove_timeout_callback(self.L["debounce_amb_factor"])
            except: pass
            
        def trigger_ui_update():
            self.L["debounce_amb_factor"] = None
            self.update_all_samples()
            
        self.L["debounce_amb_factor"] = curr_doc.add_timeout_callback(trigger_ui_update, 300)

    def on_ann_filter_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        
        is_full_range = (0 in new)
        # 1. Update UI-side state
        self.state["config"]["genome_annotations"]["show_full_range"] = is_full_range
        
        # 2. Update Backend state: 
        # Unchecked (is_full_range=False) -> Focused Extracted Range (filter_annotations=False)
        # Checked (is_full_range=True) -> Full Range (filter_annotations=True)
        if "data_selection" not in self.state["config"]: self.state["config"]["data_selection"] = {}
        target_backend_filter = is_full_range
        
        prev_filter = self.state["config"]["data_selection"].get("filter_annotations", False)
        self.state["config"]["data_selection"]["filter_annotations"] = target_backend_filter

        # If backend filtering changed, we must re-run analysis to "collect" more/less data
        if prev_filter != target_backend_filter:
            log_safe(self.state, f"Re-collecting annotations (Full Range Mode={is_full_range})...")
            with self.state["lock"]:
                self.state["analysis_data"] = None
                self.state["analysis_running"] = True
            self.L["data_rendered"] = False
            set_progress_message(self.L["div_progress"], "Re-scanning GTFs...", True)
            threading.Thread(target=run_analysis_thread, args=(self.state,), daemon=True).start()
        else:
            # Just UI refresh if backend didn't change
            self.genome_manager.on_ann_filter_change(attr, old, new)

    def on_gtf_selection_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        
        # Ensure the first GTF is always selected
        if 0 not in new:
            # Re-insert the first index and update widget (this will trigger callback again)
            self.L["sel_gtfs"].active = sorted([0] + list(new))
            return

        gtf_paths = self.state["config"]["data_selection"].get("gtf_paths", [])
        selected_paths = [gtf_paths[i] for i in new if i < len(gtf_paths)]

        self.state["config"]["data_selection"]["selected_gtfs"] = selected_paths
        log_safe(self.state, f"GTF Selection changed: {len(selected_paths)} files")
        
        # Trigger re-analysis
        with self.state["lock"]:
            self.state["analysis_data"] = None
            self.state["analysis_running"] = True
        self.L["data_rendered"] = False
        set_progress_message(self.L["div_progress"], "Re-processing GTFs...", True)
        threading.Thread(target=run_analysis_thread, args=(self.state,), daemon=True).start()

    def on_show_cohort_selection_click(self):
        gene = self.L["sel_gene"].value
        if not gene:
            return
        
        # Capture current selections to prefill the panel
        current_selection = self.L["sel_samples"].value
        
        # Update the HTML with prefilled checkboxes
        self.L["div_gene_info"].text = self._get_gene_info_html(gene, preselected_samples=current_selection)
        
        # Force visibility of the info panel even if samples are selected
        self.L["force_show_cohort_selection"] = True
        self.schedule_update_all_samples()

    def on_gene_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        
        # Explicitly clear Rust LNC1 caches and reset peak memory metrics
        import lnc_seeker
        try:
            lnc_seeker.clear_all_caches_py()
        except Exception as e:
            print(f"Warning: Failed to clear Rust caches: {e}")

        # Reset force flag when switching genes
        self.L["force_show_cohort_selection"] = False
        
        print(f"DEBUG: on_gene_change triggered with new value: {new}")
        gene = new
        
        with self.state["lock"]:
            # Reset identity mappings and clear stale data indicators
            self.state["analysis_data"] = None
            self.state["data_gene_name"] = None
            self.state["stem_to_cohort"] = {}
            self.state["cohort_to_path"] = {}
            if "data_selection" not in self.state["config"]:
                self.state["config"]["data_selection"] = {}
            self.state["config"]["data_selection"]["bam_to_cohort"] = {}

            if gene and gene in self.state.get("bam_hierarchy", {}):
                entry = self.state["bam_hierarchy"][gene]
                
                # Capture metadata range if available for explicit core range focus
                metadata = entry.get("metadata", {})
                region_str = metadata.get("region")
                if region_str and ":" in region_str and "-" in region_str:
                    try:
                        ref, coords = region_str.split(":")
                        start_s, end_s = coords.split("-")
                        self.state["config"]["data_selection"]["analysis_reference"] = ref.strip()
                        self.state["config"]["data_selection"]["analysis_start"] = int(start_s.strip())
                        self.state["config"]["data_selection"]["analysis_end"] = int(end_s.strip())
                    except Exception as e:
                        print(f"DEBUG: Failed to parse region metadata '{region_str}': {e}")

                cohorts = entry.get("cohorts", {})
                cohort_names = list(cohorts.keys())
                
                # Build identity mappings for reliable backend integration
                for c_name, c_info in cohorts.items():
                    p = c_info["path"] if isinstance(c_info, dict) else c_info
                    stem = os.path.basename(p).replace(".bam", "")
                    
                    self.state["cohort_to_path"][c_name] = p
                    self.state["config"]["data_selection"]["bam_to_cohort"][p] = c_name
                    self.state["stem_to_cohort"][stem] = c_name

                # Use shared sorting logic for consistent UI
                cohort_data = self._get_sorted_cohort_data(gene)

                sample_options = []
                for i, c in enumerate(cohort_data):
                    name, status, tissue = c["name"], c["status"], c["tissue"]
                    val = f"{i+1}. {name}"
                    label = f"{i+1}. {name} ({status}, {tissue})"
                    sample_options.append((val, label))
                
                self.L["sel_samples"].options = sample_options
                self.L["sel_samples"].value = [] # Reset selection on gene change
                
                # Sync empty selection to config
                if "data_selection" not in self.state["config"]:
                    self.state["config"]["data_selection"] = {}
                self.state["config"]["data_selection"]["selected_samples"] = []

                self.L["div_gene_info"].text = self._get_gene_info_html(gene)
            else:
                self.L["sel_samples"].options = []
                self.L["div_gene_info"].text = ""
        
        self.L["sel_samples"].value = []
        self.L["shared_rules_cache"] = {}
        self.L["src_shared_rules"].data = dict(sample=[], curated=[], predicted=[], novel=[])
        self.L["shared_rules_container"].children = []
        self.L["sample_plots"].clear()
        self.L["reset_x_range"] = True
        
        # Ensure junction selection and read caches are cleared on gene change
        self.selection_manager.on_junction_selection_change("gene_change", new=[])
        
        self.on_sample_selection_change(None, None, [])

    def on_sample_selection_change(self, attr, old, new):
        if self.state.get("config") is None:
            return
        
        # Reset the force flag whenever a selection is confirmed or changed manually
        self.L["force_show_cohort_selection"] = False
        
        self.rules_manager.update_rules_ui()
        selected_names = self.selection_manager.get_selected_samples()
        
        # Synchronize selected samples (raw names) to config for benchmarking and pipeline use
        from lnc_seeker_bokeh.state import strip_id
        self.state["config"]["data_selection"]["selected_samples"] = [strip_id(s) for s in selected_names]

        all_required = set(selected_names)
        rules_data = self.L["src_shared_rules"].data
        for r_type in ["curated", "predicted", "novel"]:
            for i, val in enumerate(rules_data.get(r_type, [])):
                if val in ["+ (Present)", "- (Absent)"]: all_required.add(rules_data['sample'][i])
        
        re_scan, re_fetch = False, False
        cur_q = self.L["sld_mq"].value
        if self.state["config"]["coverage_and_junctions_profile"].get("min_mapping_quality") != cur_q:
            self.state["config"]["coverage_and_junctions_profile"]["min_mapping_quality"] = cur_q; re_scan = True
        
        cur_gene = self.L["sel_gene"].value
        if self.state["config"]["general"].get("gene_name") != cur_gene:
            self.state["config"]["general"]["gene_name"] = cur_gene; re_scan = True
        
        cur_bams = []
        for name in all_required:
            p = self.selection_manager.get_bam_path(name)
            if p:
                cur_bams.append(p)
        
        if set(self.state["config"]["data_selection"].get("bam_paths", [])) != set(cur_bams):
            self.state["config"]["data_selection"]["bam_paths"] = cur_bams; re_scan = True
        
        if re_scan:
            should_start_thread = False
            with self.state["lock"]:
                self.state["analysis_data"] = None
                self.state["data_gene_name"] = None
                if not self.state.get("analysis_running", False):
                    self.state["analysis_running"] = True
                    should_start_thread = True
            
            self.L["data_rendered"] = False
            set_progress_message(self.L["div_progress"], "Queuing analysis...", True)
            
            if should_start_thread:
                threading.Thread(target=run_analysis_thread, args=(self.state,), daemon=True).start()
            else:
                log_safe(self.state, "Analysis already running. Coalescing request...")
        
        self.schedule_update_all_samples()

    def on_filter_flanks_change(self, attr, old, new):
        self.update_all_samples()
        self.reads_manager.update_reads_ui()

    def on_sample_multiselect_change(self, attr, old, new):
        """Debounce manual selection changes in the MultiSelect widget."""
        curr_doc = self.doc
        if self.L.get("debounce_sample_sel"):
            try:
                curr_doc.remove_timeout_callback(self.L["debounce_sample_sel"])
            except:
                pass
        
        def do_update():
            self.L["debounce_sample_sel"] = None
            self.on_sample_selection_change(attr, old, new)

        self.L["debounce_sample_sel"] = curr_doc.add_timeout_callback(do_update, 400)

    def update_normalization_ui(self):
        """Enable or disable normalization sub-options based on the master normalization toggle."""
        # Check if "Normalize by Samples" is active (id 0 in the CheckboxGroup)
        is_normalizing = (0 in self.L["chk_normalize"].active)
        self.L["chk_normalize_junctions"].disabled = not is_normalizing
        self.L["chk_global_y_range_normalize"].disabled = not is_normalizing

    def on_normalize_change(self, attr, old, new):
        """Handle toggle of 'Normalize by Samples'."""
        self.update_normalization_ui()
        self.update_all_samples()

    def _setup_callbacks(self):
        self.L["sel_gene"].on_change('value', self.benchmark.wrap_callback("Interaction", "on_gene_change")(self.on_gene_change))
        self.L["sel_samples"].on_change('value', self.benchmark.wrap_callback("Interaction", "on_sample_multiselect_change")(self.on_sample_multiselect_change))
        
        # Register range change listener for automatic resampling
        self.L["shared_x_range"].on_change('start', self.on_range_change)
        self.L["shared_x_range"].on_change('end', self.on_range_change)
        self.L["spn_height"].on_change('value', lambda a, o, n: self.update_all_samples())
        self.L["sld_mq"].on_change('value', self.benchmark.wrap_callback("Interaction", "on_mq_change")(self.on_mq_change))
        self.L["sld_min_reads"].on_change('value', lambda a, o, n: self.update_all_samples())
        self.L["chk_filter_flanks"].on_change('active', self.on_filter_flanks_change)
        self.L["chk_full_cov"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_show_bg"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_normalize"].on_change('active', self.on_normalize_change)
        self.L["chk_normalize_junctions"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_global_y_range_normalize"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_cliffs"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_ambiguity"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["sld_ambiguity_factor"].on_change('value', self.on_ambiguity_factor_change)
        self.L["spn_ann_height"].on_change('value', lambda a, o, n: self.update_all_samples())
        self.L["chk_markers"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_cds"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_unsupported_introns"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_ann_filter"].on_change('active', self.on_ann_filter_change)
        self.L["sel_gtfs"].on_change('active', self.benchmark.wrap_callback("Interaction", "on_gtf_selection_change")(self.on_gtf_selection_change))
        self.L["spn_max_reads"].on_change('value', lambda a, o, n: self.reads_manager.update_reads_ui())
        self.L["spn_reads_height"].on_change('value', lambda a, o, n: self.update_all_samples())
        self.L["sel_backend"].on_change('value', self.on_backend_change)
        self.L["chk_fixed_width"].on_change('active', self.on_fixed_width_toggle)
        self.L["spn_plot_width"].on_change('value', self.on_fixed_width_toggle)
        self.L["chk_crosshair"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_crosshair_ext"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_xaxis_main"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_legend_main"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_xaxis_ann"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_xaxis_reads"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_legend_reads"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_show_gap_size"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_show_deletions"].on_change('active', lambda a, o, n: self.update_all_samples())
        
        self.L["cohort_selection_trigger"].on_change('value', lambda a, o, n: self.on_sample_selection_change(None, None, self.L["sel_samples"].value))
        
        # Consistent registration for all buttons in the "Transcript Creator" tab
        # Using ONLY .on_click which is standard for Button widgets in Bokeh 3.x
        self.L["btn_add_j"].on_click(self.creator.on_add_j_click)
        self.L["btn_show_cohort_selection"].on_click(self.on_show_cohort_selection_click)
        self.L["btn_clear_j"].on_click(self.creator.on_clear_j_click)
        self.L["btn_remove_j"].on_click(self.creator.on_remove_j_click)
        self.L["btn_fetch_start"].on_click(self.creator.on_fetch_start_click)
        self.L["btn_fetch_end"].on_click(self.creator.on_fetch_end_click)
        self.L["btn_zoom_to_creator"].on_click(self.creator.on_zoom_to_creator_click)
        self.L["btn_create_t"].on_click(self.creator.on_create_t_click)
        
        self.L["btn_export_json"].on_click(self.creator.on_export_json_click)
        self.L["btn_export_gtf"].on_click(self.creator.on_export_gtf_click)
        self.L["btn_export_gff3"].on_click(self.creator.on_export_gff3_click)
        self.L["file_import_json"].on_change('value', self.creator.on_import_json_change)

        # Enhanced JS Trigger for multiple formats
        self.L["div_download"].js_on_change('text', CustomJS(args=dict(tid=self.L["txt_transcript_id"]), code="""
            if (!cb_obj.text || cb_obj.text.indexOf("::CONTENT::") === -1) return;
            const parts = cb_obj.text.split("::CONTENT::");
            const fmt = parts[0];
            const content = parts[1];
            
            const extensions = { 'gtf': '.gtf', 'gff3': '.gff3', 'json': '.json' };
            const mimes = { 'gtf': 'text/plain', 'gff3': 'text/plain', 'json': 'application/json' };
            
            const filename = (tid.value || "transcript") + (extensions[fmt] || ".txt");
            const blob = new Blob([content], {type: mimes[fmt] || 'text/plain'});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            cb_obj.text = ""; 
        """))

        # Bridge for selecting samples directly from the info-panel HTML via Button
        self.L["sel_gene"].js_on_change('value', CustomJS(args=dict(ms=self.L["sel_samples"], trigger=self.L["cohort_selection_trigger"]), code="""
            window.selectFromPanel = function(btn) {
                let checkboxes = [];
                const container = btn ? btn.closest('div') : null;
                if (container) {
                    checkboxes = container.querySelectorAll('.cohort-chk');
                }
                
                // Fallback 1: Search the whole document
                if (checkboxes.length === 0) {
                    checkboxes = document.querySelectorAll('.cohort-chk');
                }
                
                // Fallback 2: Deep search across Shadow DOMs (required for some Bokeh versions)
                if (checkboxes.length === 0) {
                    const findAll = (root) => {
                        let found = Array.from(root.querySelectorAll('.cohort-chk'));
                        const children = Array.from(root.querySelectorAll('*'));
                        for (const child of children) {
                            if (child.shadowRoot) {
                                found = found.concat(findAll(child.shadowRoot));
                            }
                        }
                        return found;
                    };
                    checkboxes = findAll(document.body);
                }

                const selected = [];
                checkboxes.forEach(cb => {
                    if (cb.checked) {
                        selected.push(cb.getAttribute('data-val'));
                    }
                });
                
                if (selected.length > 0) {
                    ms.value = selected;
                    // Always trigger a refresh to close the panel, even if selection didn't change
                    trigger.value = trigger.value + 1;
                } else {
                    alert("Please select at least one cohort from the list.");
                }
            };
        """))
        
        # Other event handlers
        self.L["ds_annotations"].selected.on_change('indices', self.creator.on_ann_click)
        self.L["ds_cds"].selected.on_change('indices', self.creator.on_ann_click)
        self.L["ds_transcripts"].selected.on_change('indices', self.creator.on_ann_click)
        self.L["ds_intron_markers"].selected.on_change('indices', self.creator.on_ann_click)
        self.L["mul_curr_j"].on_change('value', self.creator.on_mul_curr_j_change)

        self.L["txt_transcript_id"].on_change('value', lambda a, o, n: self.update_all_samples())
        self.L["sel_strand"].on_change('value', lambda a, o, n: self.creator.update_junction_ui())
        self.L["num_t_start"].on_change('value', lambda a, o, n: (self.creator.update_transcript_summary(), self.update_all_samples()))
        self.L["num_t_end"].on_change('value', lambda a, o, n: (self.creator.update_transcript_summary(), self.update_all_samples()))
        self.L["chk_xaxis_creator"].on_change('active', lambda a, o, n: self.update_all_samples())
        self.L["chk_markers_creator"].on_change('active', lambda a, o, n: self.update_all_samples())
        
        self.L["chk_types"].on_change('active', self.on_junction_types_change)
        self.L["tabs"].on_change('active', lambda a, o, n: self.update_all_samples())

        curr_doc = self.doc
        curr_doc.add_periodic_callback(self.update_progress, 500)
        
        self.doc.add_root(self.layout)
