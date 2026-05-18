# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import math
import json
import threading
import time
from bokeh.plotting import curdoc
from bokeh.models import (
    Range1d, PanTool, WheelZoomTool, CustomJS, 
    NumeralTickFormatter, LabelSet, CustomAction
)
from bokeh.events import Reset, MouseWheel
import lnc_seeker
from lnc_seeker_bokeh.state import log_safe
from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot
from lnc_seeker_bokeh.constants import (
    set_progress_message, clear_progress, 
    set_progress_complete, set_progress_fail,
    set_progress_success
)
from lnc_seeker_bokeh.reads_plot import create_reads_plot, apply_reads_visibility, apply_reads_styles

class ReadsManager:
    def __init__(self, app):
        """
        Coordinates the 'Full Read Layout' plot and data fetching.
        :param app: The VisualizerApp instance this manager belongs to.
        """
        self.app = app
        self.L = app.L
        self.request_id = 0

    def _get_fetch_key(self, selected_junctions, target_samples_bams):
        """Generates a unique key for the current fetch parameters to detect redundant requests."""
        L = self.L
        state = self.app.state
        
        # Core parameters that determine BAM fetch results
        mq = state["config"]["coverage_and_junctions_profile"].get("min_mapping_quality", 20)
        max_reads = L["spn_max_reads"].value if "spn_max_reads" in L else 100
        filter_clean = (0 in L["chk_filter_flanks"].active) if "chk_filter_flanks" in L else False
        genome_path = state["config"]["data_selection"].get("genome_path", "")
        
        # Return a stable tuple of all relevant factors
        return (
            tuple(sorted(selected_junctions)), 
            tuple(sorted(target_samples_bams)),
            mq, 
            max_reads, 
            filter_clean, 
            genome_path
        )

    def update_reads_ui(self, selected_junctions=None, target_samples_bams=None, cached_data=None):
        """
        Unified entry point for update full-read layouts. 
        Always uses the 'collective' logic (even for 1 sample/1 junction).
        """
        # If no arguments provided, try to restore from last known state
        if selected_junctions is None:
            selected_junctions = self.L.get("last_selected_junctions", [])
        if target_samples_bams is None:
            target_samples_bams = self.L.get("last_target_samples_bams", [])
            
        if not selected_junctions or not target_samples_bams:
            if cached_data is None: 
                return

        # Core logic cleanup: remove obsolete state variables
        self.L["last_selected_junctions"] = selected_junctions
        self.L["last_target_samples_bams"] = target_samples_bams
        
        # Check if we can reuse the existing collective payload cache
        is_major_change = True
        if cached_data is None:
            current_key = self._get_fetch_key(selected_junctions, target_samples_bams)
            last_key = self.L.get("last_read_fetch_key")
            
            if last_key is not None:
                # Key structure: (junctions, samples, mq, max_reads, filter_clean, genome_path)
                # Major changes: junctions or samples
                is_major_change = current_key[:2] != last_key[:2]

            if current_key == last_key and self.L.get("cached_collective_payload") is not None:
                # Parameters haven't changed, reuse the loaded data WITHOUT refetching from BAM
                cached_data = self.L["cached_collective_payload"]
            else:
                # Update the key for the new upcoming fetch
                self.L["last_read_fetch_key"] = current_key

        doc = self.L.get("doc") or curdoc()
        
        # Set busy flag to drop scroll events during drawing
        if "ds_reads_busy" in self.L:
            self.L["ds_reads_busy"].data = dict(is_busy=[1])

        # Magic-key logic: version individual fetch requests
        self.request_id += 1
        curr_id = self.request_id

        # Set immediate notification on main thread only if actually fetching
        if cached_data is None:
            self.L["is_fetching_reads"] = True
            num_j = len(selected_junctions)
            num_s = len(target_samples_bams)
            msg = f"Fetching reads for {num_j} junction{'s' if num_j > 1 else ''} in {num_s} sample{'s' if num_s > 1 else ''}..."
            set_progress_message(self.L["div_progress"], msg, True)
        
        def run_fetch():
            try:
                # Small sleep ensures the UI message is flushed before GIL-heavy work
                time.sleep(0.01)
                self._fetch_and_render_reads(selected_junctions, target_samples_bams, cached_data, curr_id, is_major_change)
            except Exception as e:
                import traceback
                log_safe(self.app.state, f"Error fetching reads: {e}\n{traceback.format_exc()}")
                if curr_id == self.request_id:
                    self.L["is_fetching_reads"] = False
                    doc.add_next_tick_callback(lambda: set_progress_fail(self.L["div_progress"], f"Error fetching reads: {str(e)[:50]}..."))
                    # Fallback reset of the busy flag if drawing fails
                    if "ds_reads_busy" in self.L:
                        doc.add_next_tick_callback(lambda: self.L["ds_reads_busy"].data.update(dict(is_busy=[0])))
        
        def start_fetch():
            if cached_data is not None:
                run_fetch()
            else:
                threading.Thread(target=run_fetch, daemon=True).start()

        # Push to next tick to ensure set_progress_message is registered by the browser
        doc.add_next_tick_callback(start_fetch)

    def _fetch_and_render_reads(self, selected_junctions, target_samples_bams, cached_data=None, req_id=None, is_major_change=True):
        """Background logic for fetching and calculating read positions."""
        L = self.L
        state = self.app.state
        doc = L.get("doc") or curdoc()
        
        all_reads_count = 0
        ex_x0, ex_x1, ex_y, ex_colors, ex_names, ex_mismatches, ex_insertions, ex_types, ex_thickness, ex_labels, ex_mid_x, ex_label_colors, ex_label_bg_colors = [], [], [], [], [], [], [], [], [], [], [], [], []
        in_x0, in_x1, in_y, in_labels, in_mid_x, in_colors, in_widths, in_names, in_types, in_mismatches, in_insertions = [], [], [], [], [], [], [], [], [], [], []
        br_x0, br_x1, br_y, br_labels, br_mid_x, br_names = [], [], [], [], [], []
        full_x0, full_x1, full_y, full_names = [], [], [], []
        section_x, section_y, section_text = [], [], []

        # Use fixed step for vertical layout, letting Bokeh handle zooming/squeezing via y_range
        v_step = 1.0
        r_styles = state["config"]["full_read_layout"].get("layout_styles", {})
        filter_clean = (0 in L["chk_filter_flanks"].active) if "chk_filter_flanks" in L else False
        
        if cached_data is not None:
            all_sample_results = cached_data
        else:
            t_fetch_start = time.time()
            all_sample_results = []
            
            # Group junctions by reference for batch fetching
            j_by_ref = {}
            for j in selected_junctions:
                ref = str(j[0])
                if ref not in j_by_ref: j_by_ref[ref] = []
                j_by_ref[ref].append((int(j[1]), int(j[2])))

            mq = state["config"]["coverage_and_junctions_profile"].get("min_mapping_quality", 20)
            max_r = int(self.L["spn_max_reads"].value if "spn_max_reads" in self.L else 100)
            g_path = state["config"]["data_selection"].get("genome_path")

            for s_name, t_bam in target_samples_bams:
                sample_reads = []
                seen_read_names = set()
                
                # Fetch reads for each reference chromosome involved
                for ref_name, j_tuples in j_by_ref.items():
                    try:
                        reads_batch_json = lnc_seeker.get_junction_reads_batch_py(
                            t_bam, ref_name, j_tuples, mq, max_r, g_path, filter_clean
                        )
                        batch_data = json.loads(reads_batch_json)
                        for _j_key, reads in batch_data.items():
                            for r in reads:
                                if r['name'] not in seen_read_names:
                                    sample_reads.append(r)
                                    seen_read_names.add(r['name'])
                    except Exception as e:
                        log_safe(state, f"Batch fetch failed for {ref_name} in {s_name}: {e}")

                if sample_reads:
                    all_sample_results.append((s_name, sample_reads))
            
            fetch_duration = time.time() - t_fetch_start
            fetch_details = f"{len(selected_junctions)} juncs, {len(target_samples_bams)} samples"
            
            # Use explicit probe methods to fill metrics
            c_core = self.app.get_cache_core_mb()
            c_related = self.app.get_cache_related_mb()
            c_annotation = self.app.get_cache_annotation_mb()
            c_json = self.app.get_json_overhead_mb()

            self.app.benchmark.log_event("Interaction", "collect_supporting_reads", fetch_duration, 
                                       reads=sum(len(r) for _, r in all_sample_results),
                                       cache_core_mb=c_core, cache_related_mb=c_related,
                                       cache_annotation_mb=c_annotation, json_overhead_mb=c_json,
                                       details=fetch_details)
            
            # Prioritize reads that actually span any of the target junctions (Green reads to top)
            target_juncs = set((int(j[1]), int(j[2])) for j in selected_junctions)
            def is_target_read(r):
                m_map = {}
                for seg in r.get('segments', []):
                    m_id = seg.get('is_mate', False)
                    if m_id not in m_map: m_map[m_id] = []
                    m_map[m_id].append(seg)
                for m_id, segs in m_map.items():
                    segs = sorted(segs, key=lambda x: x['start'])
                    for i in range(len(segs)-1):
                        if (int(segs[i]['end']), int(segs[i+1]['start'])) in target_juncs: return True
                return False

            for i in range(len(all_sample_results)):
                s_name, s_reads = all_sample_results[i]
                # Sort by: 1. Spans target junction, 2. Extension size (end - start)
                all_sample_results[i] = (s_name, sorted(s_reads, key=lambda r: (is_target_read(r), int(r['end']) - int(r['start'])), reverse=True))
            
            L["cached_collective_payload"] = all_sample_results

        total_h = 0
        for s_name, reads in all_sample_results:
            n = len(reads)
            all_reads_count += n
            total_h += n * v_step + (v_step * 7.0)

        if not selected_junctions: 
            # Reset busy flag if we abort early
            if "ds_reads_busy" in self.L:
                doc.add_next_tick_callback(lambda: self.L["ds_reads_busy"].data.update(dict(is_busy=[0])))
            return 
            
        ref_main = selected_junctions[0][0]
        start_min = min(j[1] for j in selected_junctions)

        current_y = total_h
        for s_name, reads in all_sample_results:
            current_y -= (v_step * 5.0)
            section_x.append(start_min)
            section_y.append(current_y)
            
            # Fetch metadata for the section label
            tissue, status, num_samples = self.app.selection_manager.get_sample_metadata(s_name)
            section_text.append(f"{s_name} ({tissue}, {status}, samples={num_samples})")
            
            current_y -= v_step * 2.0 # Balanced gap

            for read in reads:
                y_val = current_y
                full_x0.append(read['start'])
                full_x1.append(read['end'])
                full_y.append(y_val)
                full_names.append(read['name'])
                
                mate_map = {}
                for seg in read['segments']:
                    m_id = seg.get('is_mate', False)
                    if m_id not in mate_map: mate_map[m_id] = []
                    mate_map[m_id].append(seg)
                
                introns_by_mate = {}
                for m_id, m_segments in mate_map.items():
                    m_segments = sorted(m_segments, key=lambda x: x['start'])
                    introns_by_mate[m_id] = []
                    for k, seg in enumerate(m_segments):
                        ex_x0.append(seg['start']); ex_x1.append(seg['end']); ex_y.append(y_val); ex_names.append(read['name'])
                        mism = seg.get('mismatches', 0)
                        ins = seg.get('insertions', 0)

                        ex_mismatches.append(mism)
                        ex_insertions.append(ins)
                        ex_types.append("exon")
                        ex_mid_x.append((seg['start'] + seg['end']) / 2.0)
                        
                        label_parts = []
                        if mism > 0: label_parts.append(f"X{mism}")
                        if ins > 0: label_parts.append(f"I{ins}")
                        ex_labels.append(" ".join(label_parts))
                        
                        # Use base styles and white font in colored boxes
                        base_lw = r_styles.get("show_reads_line_width", 2)
                        style_cfg = r_styles.get("exon_mate" if seg.get("is_mate") else "exon_single", {})
                        color = style_cfg.get("color", "royalblue")
                        if (mism > 0 or ins > 0) and "variant_color" in style_cfg:
                            color = style_cfg["variant_color"]
                        
                        ex_colors.append(color)
                        ex_thickness.append(base_lw)
                        ex_label_colors.append("white")
                        ex_label_bg_colors.append("red" if (mism > 2 or ins > 2) else "black")

                        if k < len(m_segments) - 1:
                            next_seg = m_segments[k+1]
                            is_del = seg.get('is_followed_by_deletion', False)
                            is_clean = (mism == 0 and ins == 0 and 
                                        next_seg.get('mismatches', 0) == 0 and 
                                        next_seg.get('insertions', 0) == 0)
                            introns_by_mate[m_id].append((int(seg['end']), int(next_seg['start']), is_del, is_clean))

                all_pair_introns = {}
                for m_id, introns in introns_by_mate.items():
                    for j_start, j_end, is_del, is_clean in introns:
                        if filter_clean and not is_clean: continue
                        key = (j_start, j_end, is_del)
                        all_pair_introns[key] = all_pair_introns.get(key, 0) + 1
                
                for (j_start, j_end, is_del), count in all_pair_introns.items():
                    in_x0.append(j_start); in_x1.append(j_end); in_y.append(y_val); in_labels.append("2x" if count > 1 else ""); in_mid_x.append((j_start + j_end) / 2.0); in_names.append(read['name'])
                    in_types.append("deletion" if is_del else "junction")
                    in_mismatches.append(0)
                    in_insertions.append(0)
                    if is_del:
                        in_colors.append(r_styles.get("deletion", {}).get("color", "red"))
                        in_widths.append(r_styles.get("deletion", {}).get("line_width", 2))
                    else:
                        is_p_target = any(j_start == int(j[1]) and j_end == int(j[2]) for j in selected_junctions)
                        in_colors.append(r_styles.get("intron_target" if is_p_target else "intron_other", {}).get("color", "green"))
                        in_widths.append(r_styles.get("intron_target" if is_p_target else "intron_other", {}).get("line_width", 2 if is_p_target else 1))

                # Connect mates with bridges
                if len(mate_map) > 1:
                    m_ids = sorted(mate_map.keys(), key=lambda m: min(s['start'] for s in mate_map[m]))
                    for i in range(len(m_ids) - 1):
                        m1_id, m2_id = m_ids[i], m_ids[i+1]
                        m1_end = max(s['end'] for s in mate_map[m1_id])
                        m2_start = min(s['start'] for s in mate_map[m2_id])
                        if m2_start > m1_end:
                            br_x0.append(m1_end); br_x1.append(m2_start); br_y.append(y_val); br_names.append(read['name'])
                            br_labels.append(f"mate gap: {m2_start - m1_end}bp"); br_mid_x.append((m1_end + m2_start) / 2.0)

                current_y -= v_step
        
        def commit_to_bokeh():
            if req_id is not None and req_id != self.request_id: return
            L["ds_reads_exons"].data = dict(x0=ex_x0, x1=ex_x1, y=ex_y, color=ex_colors, thickness=ex_thickness, 
                                            name=ex_names, mismatches=ex_mismatches, insertions=ex_insertions, 
                                            type=ex_types, label=ex_labels, mid_x=ex_mid_x, label_color=ex_label_colors,
                                            label_bg_color=ex_label_bg_colors)
            L["ds_reads_introns"].data = dict(x0=in_x0, x1=in_x1, y=in_y, label=in_labels, mid_x=in_mid_x, color=in_colors, line_w=in_widths, name=in_names, type=in_types, mismatches=in_mismatches, insertions=in_insertions)
            L["ds_reads_bridges"].data = dict(x0=br_x0, x1=br_x1, y=br_y, label=br_labels, mid_x=br_mid_x, name=br_names, color=[""] * len(br_x0), type=["bridge"]*len(br_x0), mismatches=[0]*len(br_x0), insertions=[0]*len(br_x0))
            L["ds_reads_full"].data = dict(x0=full_x0, x1=full_x1, y=full_y, name=full_names, color=["#444444"] * len(full_x0))
            L["ds_reads_labels"].data = dict(x=section_x, y=section_y, text=section_text, color=["white"] * len(section_x))
            
            if L["p_reads"]:
                ranges_str = ", ".join([f"{int(j[1])}-{int(j[2])}" for j in selected_junctions[:3]])
                L["p_reads"].title.text = f"Full Read Layouts ({all_reads_count}) - {ref_main}:{ranges_str}"
                
                # Using standard Zoom with Range1d, similar to Genome Annotations.
                # Only reset range on full selection changes (not cached redraws).
                if cached_data is None:
                    last_h = L.get("last_reads_total_h")
                    # Preserving the squeezing-level (span) across both major and minor changes
                    current_span = L["p_reads"].y_range.end - L["p_reads"].y_range.start
                    
                    if is_major_change or last_h is None:
                        # Scroll to top for major changes, but keep current zoom level
                        L["p_reads"].y_range.start = total_h + v_step - current_span
                        L["p_reads"].y_range.end = total_h + v_step
                    else:
                        # Shift existing range to maintain position when reads are filtered
                        delta_h = total_h - last_h
                        L["p_reads"].y_range.start += delta_h
                        L["p_reads"].y_range.end += delta_h
                
                L["last_reads_total_h"] = total_h
            
            if cached_data is None:
                L["is_fetching_reads"] = False
                L["is_sticky_message"] = True
                set_progress_success(L["div_progress"], f"Layout updated with {all_reads_count} reads.")
                def clear(): 
                    L["is_sticky_message"] = False
                    doc.add_next_tick_callback(lambda: clear_progress(L["div_progress"]))
                threading.Timer(2.0, clear).start()
            
            # Reset busy flag (legacy support for drawing events)
            if "ds_reads_busy" in self.L:
                self.L["ds_reads_busy"].data = dict(is_busy=[0])
        doc.add_next_tick_callback(commit_to_bokeh)

    def get_or_create_reads_plot(self, reads_target_h, fixed_width, target_width, backend):
        """Ensures the reads plot exists and is configured correctly."""
        if self.L.get("p_reads") is None:
            p_reads = create_reads_plot(
                self.L, self.app.state, reads_target_h, 
                fixed_width, target_width, backend, self.app.on_reset_zoom_click
            )
            self.L["p_reads"] = p_reads
        else:
            p_reads = self.L["p_reads"]
            if p_reads.x_range != self.L["shared_x_range"]:
                p_reads.x_range = self.L["shared_x_range"]
            p_reads.height = reads_target_h
            apply_reads_visibility(p_reads, self.L, self.app.state)
            p_reads.output_backend = backend

            if fixed_width:
                p_reads.sizing_mode = "fixed"; p_reads.width = target_width
            else:
                p_reads.sizing_mode = "stretch_width"
            
            r_styles_ru = self.app.state["config"]["full_read_layout"].get("layout_styles", {})
            if r_styles_ru:
                apply_reads_styles(self.L, r_styles_ru)
        
        return p_reads
