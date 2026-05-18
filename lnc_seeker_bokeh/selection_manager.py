# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import threading
from bokeh.plotting import curdoc
from lnc_seeker_bokeh.state import log_safe, strip_id
from lnc_seeker_bokeh.constants import set_progress_message, clear_progress, set_progress_fail

class SelectionManager:
    def __init__(self, app):
        self.app = app
        self.L = app.L

    def get_raw_name(self, name):
        """Standardized utility to get the clean cohort name without UI prefix."""
        return strip_id(name)

    def get_selected_samples(self):
        """Returns the currently selected samples sorted by their order in the MultiSelect options."""
        selected = self.L["sel_samples"].value
        options = self.L["sel_samples"].options
        
        if not options:
            return list(selected)
            
        # Create a mapping of value to index from options
        opt_values = {}
        for i, opt in enumerate(options):
            if isinstance(opt, (list, tuple)):
                opt_values[opt[0]] = i
            else:
                opt_values[opt] = i
        
        # Sort based on the index in options
        return sorted(list(selected), key=lambda x: opt_values.get(x, 999))

    def get_sample_metadata(self, name):
        """Retrieves tissue, status, and num_samples for a given sample/cohort name."""
        raw_name = strip_id(name)
        gene = self.L["sel_gene"].value if "sel_gene" in self.L else None
        
        tissue, status, num_samples = "N/A", "N/A", "N/A"
        
        if gene and gene in self.app.state.get("bam_hierarchy", {}):
            cohorts = self.app.state["bam_hierarchy"][gene].get("cohorts", {})
            if raw_name in cohorts:
                c_info = cohorts[raw_name]
                if isinstance(c_info, dict):
                    tissue = c_info.get("tissue", "N/A")
                    status = c_info.get("status", "N/A")
                    num_samples = str(c_info.get("num_samples", "N/A"))
        
        return tissue, status, num_samples

    def get_bam_path(self, name):
        """Maps a UI name or raw cohort name to its absolute BAM path."""
        raw = strip_id(name)
        return self.app.state.get("cohort_to_path", {}).get(raw)

    def get_sample_data(self, name):
        """Retrieves results for a sample from the current analysis state."""
        raw = strip_id(name)
        data = self.app.state.get("analysis_data", {})
        if data:
            return data.get("samples", {}).get(raw)
        return None

    def on_intron_selection_change(self, attr, old, new):
        """Handles selecting an intron in the genome annotations plot."""
        if self.L.get("selection_updating", False):
            return
        
        ds_introns = self.L["ds_introns"]
        
        # Collect all selected introns as (ref, start, end) tuples to trigger 
        # synchronization across all components (plots, creator preview, etc.)
        all_selected_tuples = []
        found_experimental_match = False
        
        if new:
            for idx in new:
                target_start = int(ds_introns.data['left'][idx])
                target_end = int(ds_introns.data['right'][idx])
                target_ref = str(ds_introns.data['reference'][idx])
                all_selected_tuples.append((target_ref, target_start, target_end))
                
                # Check for experimental matches just for the warning notification
                if not found_experimental_match:
                    for s_name, plot_objs in self.L["sample_plots"].items():
                        ds_js = plot_objs[2]
                        if not ds_js.data or 'x0' not in ds_js.data:
                            continue
                        # Search for exact coordinate match in experimental data
                        d_js = ds_js.data
                        for i in range(len(d_js['x0'])):
                            if (str(d_js['reference'][i]) == target_ref and 
                                int(d_js['x0'][i]) == target_start and 
                                int(d_js['x1'][i]) == target_end):
                                found_experimental_match = True
                                break
                        if found_experimental_match: 
                            break

        # Always trigger master sync - this will now also update the Transcript Creator
        self.on_junction_selection_change("intron", new=all_selected_tuples)
        
        if new and not found_experimental_match:
            # No experimental match found in any diagram - show warning but keep the selection
            log_safe(self.app.state, "No exact experimental match found for selected intron.")
            
            # Visual notification
            self.L["is_sticky_message"] = True
            set_progress_fail(self.L["div_progress"], "No experimental junction matching intron found.")
            
            def delayed_clear():
                self.L["is_sticky_message"] = False
                doc = self.L.get("doc") or curdoc()
                doc.add_next_tick_callback(lambda: clear_progress(self.L["div_progress"]))
            
            threading.Timer(3.0, delayed_clear).start()

    def clear_read_layout(self):
        """Specifically clears the read layout sources."""
        empty_data = dict(x0=[], x1=[], y=[], color=[], name=[])
        empty_intron = dict(x0=[], x1=[], y=[], label=[], mid_x=[], color=[], line_w=[], name=[])
        empty_bridge = dict(x0=[], x1=[], y=[], label=[], mid_x=[], name=[], color=[])
        
        self.L["ds_reads_exons"].data = empty_data
        self.L["ds_reads_introns"].data = empty_intron
        self.L["ds_reads_bridges"].data = empty_bridge
        self.L["ds_reads_full"].data = empty_data
        self.L["ds_reads_labels"].data = dict(x=[], y=[], text=[], color=[])

        # Clear read cache to ensure fresh fetch when new junctions are clicked
        self.L["cached_collective_payload"] = None
        self.L["last_read_fetch_key"] = None
        
        if self.L.get("p_reads"):
            self.L["p_reads"].title.text = "Click a junction label to see full read layout"

    def on_junction_selection_change(self, source, attr=None, old=None, new=None):
        """Handles selecting one or more junctions in any of the sample plots."""
        if self.L.get("selection_updating", False): 
            return
        
        self.L["selection_updating"] = True
        
        try:
            selected_tuples = []
            trigger_indices = []
            found_trigger_plot = None
            
            # 1. Identify what was selected
            if source in ["sample_plot", "creator", "intron"]:
                selected_tuples = new
            else:
                # Sync selection across all plots
                for s_name, plot_objs in self.L["sample_plots"].items():
                    ds_js = plot_objs[2]
                    if ds_js == source:
                        for idx in new:
                            d = ds_js.data
                            selected_tuples.append((str(d['reference'][idx]), int(d['x0'][idx]), int(d['x1'][idx])))
                        found_trigger_plot = ds_js
                        trigger_indices = new
                        break
                    
            if not selected_tuples:
                self.L["last_selected_junctions"] = []
                log_safe(self.app.state, "Selection cleared.")
                self.clear_read_layout()
            else:
                self.L["last_selected_junctions"] = selected_tuples
                log_safe(self.app.state, f"Stored {len(selected_tuples)} junctions in last_selected_junctions.")

            # 2. Update visuals for selected junctions in all plots
            for s_name, plot_objs in self.L["sample_plots"].items():
                ds_js = plot_objs[2]
                if ds_js.data:
                    data = ds_js.data
                    new_fs, new_fst, new_lw, new_lwid, new_lhei = [], [], [], [], []
                    indices = []
                    base_lw = 2
                    
                    for i in range(len(data.get('x0', []))):
                        try:
                            reads = float(data['label'][i])
                        except (ValueError, TypeError):
                            reads = 0.0
                        ref = str(data['reference'][i])
                        start, end = int(data['x0'][i]), int(data['x1'][i])
                        
                        if (ref, start, end) in selected_tuples:
                            indices.append(i)
                            new_lwid.append(30 * (2 if reads > 99 else 1))
                            new_lhei.append(20)
                            new_fs.append("12pt")
                            new_fst.append("bold")
                            new_lw.append(base_lw + 5)
                            if found_trigger_plot is None:
                                found_trigger_plot = ds_js
                                trigger_indices = indices
                        else:
                            new_lwid.append(20 * (2 if reads > 99 else 1))
                            new_lhei.append(14)
                            new_fs.append("8pt")
                            new_fst.append("normal")
                            new_lw.append(base_lw)
                    
                    ds_js.data.update({
                        'fsize': new_fs, 'fstyle': new_fst, 'line_w': new_lw, 
                        'label_w': new_lwid, 'label_h': new_lhei
                    })
                    ds_js.trigger('data', ds_js.data, ds_js.data)
                    
                    # Ensure indices are also synced across all plots to prevent accidental clearing
                    # from plots that don't have the current selection active in Bokeh.
                    if list(ds_js.selected.indices) != indices:
                        ds_js.selected.indices = indices
                
                if not selected_tuples and not found_trigger_plot:
                    found_trigger_plot = ds_js
                    trigger_indices = []

            # 3. Synchronize Bokeh selection state if needed
            # Use "sample_plot" as the source indicator for internal sync calls
            if found_trigger_plot is not None and source not in ["sample_plot", "intron", "creator"]:
                self.L["selection_updating"] = False
                # Re-trigger with a generic string source to signal this is a sync operation
                self.on_junction_selection_change("sample_plot", new=selected_tuples)
                self.L["selection_updating"] = True
            
            # 4. Synchronize Genome Annotations selection
            ds_introns = self.L.get("ds_introns")
            if ds_introns and ds_introns.data and 'left' in ds_introns.data:
                ann_indices = []
                d_ann = ds_introns.data
                for i in range(len(d_ann['left'])):
                    ref = d_ann.get('reference', [None] * len(d_ann['left']))[i]
                    if (ref, int(d_ann['left'][i]), int(d_ann['right'][i])) in selected_tuples:
                        ann_indices.append(i)
                
                if list(ds_introns.selected.indices) != ann_indices:
                    ds_introns.selected.indices = ann_indices

            # 5. Fetch reads for the selection - only if this was an actual UI interaction trigger
            # or from component logic (not a secondary sync call)
            if selected_tuples and source != "sample_plot":
                target_samples_bams = []
                for s_name in self.get_selected_samples():
                    p = self.get_bam_path(s_name)
                    if p:
                        target_samples_bams.append((s_name, p))
                            
                self.app.reads_manager.update_reads_ui(selected_tuples, target_samples_bams)

            # 6. Synchronize Transcript Creator selection
            ds_creator = self.L.get("ds_creator_annotations")
            if ds_creator and ds_creator.data and 'left' in ds_creator.data:
                creator_indices = []
                d_creator = ds_creator.data
                mul_select_vals = []
                # We match based on (start, end) coordinates
                selected_coords = set((t[1], t[2]) for t in selected_tuples)
                
                for i in range(len(d_creator['left'])):
                    if d_creator['feature'][i] == 'Intron':
                        js, je = int(d_creator['left'][i]), int(d_creator['right'][i])
                        if (js, je) in selected_coords:
                            creator_indices.append(i)
                            mul_select_vals.append(f"{js}-{je}")
                
                # Highlight in preview plot
                if list(ds_creator.selected.indices) != creator_indices:
                    ds_creator.selected.indices = creator_indices
                
                # Update sidebar MultiSelect list
                if "mul_curr_j" in self.L and self.L["mul_curr_j"]:
                    if set(self.L["mul_curr_j"].value) != set(mul_select_vals):
                         self.L["mul_curr_j"].value = mul_select_vals
                
                # Update info div for the creator
                if creator_indices and "div_creator_info" in self.L:
                    idx = creator_indices[0]
                    text = f"<b>Selection:</b> {d_creator['feature'][idx]} {d_creator['feature_num'][idx]} | "
                    text += f"Location: {int(d_creator['left'][idx]):,} - {int(d_creator['right'][idx]):,} | "
                    text += f"Size: {int(d_creator['width'][idx]):,} bp"
                    self.L["div_creator_info"].visible = True
                    self.L["div_creator_info"].text = text
                elif "div_creator_info" in self.L:
                    self.L["div_creator_info"].visible = False
        finally:
            self.L["selection_updating"] = False
