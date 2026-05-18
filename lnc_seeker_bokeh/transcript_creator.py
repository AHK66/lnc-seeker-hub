# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import sys
import lnc_seeker
from lnc_seeker_bokeh.state import log_safe
from lnc_seeker_bokeh.data_utils import process_analysis_data

class TranscriptCreator:
    def __init__(self, app):
        """
        Coordinates the 'Transcript Creator' tab logic.
        :param app: The VisualizerApp instance this creator belongs to.
        """
        self.app = app
        self.L = app.L
        self.session_id = app.session_id

    def update_junction_ui(self):
        """Updates the sidebar list of junctions and the summary div."""
        strand = self.L["sel_strand"].value if "sel_strand" in self.L else "+"
        # Sort based on strand: ascending for +, descending for -
        j_list = sorted(self.L["user_junctions"], reverse=(strand == "-"))
        
        # Use (value, label) tuples for MultiSelect options. 
        # Value remains 'start-end' for easy parsing, label adds the index.
        opts = [(f"{j[0]}-{j[1]}", f"{i+1}: {j[0]}-{j[1]}") for i, j in enumerate(j_list)]
        self.L["mul_curr_j"].options = opts
        self.L["mul_curr_j"].title = f"Current Junctions: ({len(opts)})"
        self.update_transcript_summary()
        self.app.update_all_samples()

    def update_transcript_summary(self):
        """Calculates and displays a summary of the transcript being built."""
        strand = self.L["sel_strand"].value if "sel_strand" in self.L else "+"
        junctions = sorted(self.L["user_junctions"], reverse=(strand == "-"))
        
        box_style = (
            "background-color: #f8fbff; "
            "border: 1px solid #d0e7ff; "
            "border-radius: 6px; "
            "padding: 12px; "
            "margin-top: 10px; "
            "font-family: inherit; "
            "color: #2c3e50; "
            "line-height: 1.4;"
        )

        if not junctions:
            self.L["div_transcript_summary"].text = f"""
            <div style="{box_style}">
                <b style="color: #339af0;">Transcript Summary</b><br>
                <span style="color: #868e96; font-style: italic; font-size: 0.9em;">No junctions added.</span>
            </div>
            """
            return
        
        s = self.L["num_t_start"].value
        e = self.L["num_t_end"].value
        
        num_exons = len(junctions) + 1
        total_len = 0
        genomic_range_str = "Incomplete (Set Start/End)"
        
        if s is not None and e is not None:
            # Length calculation
            # Exon 1: s to junctions[0][0] (if +) or junctions[0][1] to e (if -)
            if strand == "+":
                total_len = (junctions[0][0] - s)
            else:
                total_len = (e - junctions[0][1])
                
            # Internal Exons
            for i in range(len(junctions)-1):
                if strand == "+":
                    total_len += (junctions[i+1][0] - junctions[i][1])
                else:
                    total_len += (junctions[i][0] - junctions[i+1][1])
            
            # Last Exon: junctions[-1][1] to e (if +) or s to junctions[-1][0] (if -)
            if strand == "+":
                total_len += (e - junctions[-1][1])
            else:
                total_len += (junctions[-1][0] - s)
            
            genomic_range_str = f"{min(s, e):,} - {max(s, e):,}"
        
        # Build Overview Section
        summary_html = f"""
        <div style="{box_style}">
            <div style="margin-bottom: 8px; border-bottom: 1px solid #d0e7ff; padding-bottom: 5px;">
                <b style="font-size: 1.1em; color: #339af0;">Transcript Overview</b>
            </div>
            <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                <tr><td><b>Strand:</b></td><td style="text-align: right;">{strand}</td></tr>
                <tr><td><b>Exons:</b></td><td style="text-align: right;">{num_exons}</td></tr>
                <tr><td><b>Total Length:</b></td><td style="text-align: right;"><b>{total_len:,} bp</b></td></tr>
            </table>
            <div style="margin-top: 5px; font-size: 0.85em; color: #495057;">
                <b>Genomic:</b> {genomic_range_str}
            </div>

            <details style="margin-top: 10px; border-top: 1px dashed #ced4da; padding-top: 8px;">
                <summary style="cursor: pointer; color: #444; font-weight: bold; font-size: 0.85em;">
                    Structure Details (Exons/Introns)
                </summary>
                <div style="margin-top: 8px; font-size: 0.85em; max-height: 180px; overflow-y: auto; padding-right: 5px; border-left: 2px solid #e7f5ff; padding-left: 8px;">
        """

        if strand == "+":
            if s is not None: summary_html += f"<b>Start (5'):</b> {s:,}<br>"
            for i, (js, je) in enumerate(junctions):
                summary_html += f"<span style='color: #495057;'>Exon {i+1} End:</span> {js:,}<br>"
                summary_html += f"<div style='color: #adb5bd; margin: 2px 0; border-left: 1px dashed #adb5bd; padding-left: 5px; font-style: italic;'>Intron {i+1}</div>"
                summary_html += f"<span style='color: #495057;'>Exon {i+2} Start:</span> {je:,}<br>"
            if e is not None: summary_html += f"<b>End (3'):</b> {e:,}<br>"
        else:
            if e is not None: summary_html += f"<b>Start (5'):</b> {e:,}<br>"
            for i, (js, je) in enumerate(junctions):
                summary_html += f"<span style='color: #495057;'>Exon {i+1} Start:</span> {je:,}<br>"
                summary_html += f"<div style='color: #adb5bd; margin: 2px 0; border-left: 1px dashed #adb5bd; padding-left: 5px; font-style: italic;'>Intron {i+1}</div>"
                summary_html += f"<span style='color: #495057;'>Exon {i+2} End:</span> {js:,}<br>"
            if s is not None: summary_html += f"<b>End (3'):</b> {s:,}<br>"

        summary_html += """
                </div>
            </details>
        </div>
        """
        
        self.L["div_transcript_summary"].text = summary_html

    def on_add_j_click(self, event=None):
        """Callback for 'Add Selected Junctions' button."""
        try:
            log_safe(self.app.state, "Python: Adding selected junctions...")
            
            
            restored = self.L.get("last_selected_junctions", [])
            if not restored:
                log_safe(self.app.state, "Python: No junctions selected. Click junction labels in plots first.")
                return

            # Logic to remove all overlapping introns in the list before adding the new ones.
            # This guarantees a clean intron layout.
            current_junctions = self.L["user_junctions"]
            new_selection = [(s, e) for (ref, s, e) in restored]
            
            # Identify junctions to keep (those that don't overlap with ANY of the new ones)
            to_keep = []
            removed_count = 0
            for (cs, ce) in current_junctions:
                overlap = False
                for (ns, ne) in new_selection:
                    # Standard overlap check: (start1 < end2) and (end1 > start2)
                    if cs < ne and ce > ns:
                        overlap = True
                        break
                if not overlap:
                    to_keep.append((cs, ce))
                else:
                    removed_count += 1
            
            if removed_count > 0:
                log_safe(self.app.state, f"Python: Removed {removed_count} overlapping junctions.")

            # Add the new ones (avoid duplicates just in case)
            added = 0
            for (ns, ne) in new_selection:
                if (ns, ne) not in to_keep:
                    to_keep.append((ns, ne))
                    added += 1
            
            self.L["user_junctions"] = to_keep
            
            if added > 0:
                log_safe(self.app.state, f"Python: Added {added} junctions.")
                self.update_junction_ui()
            else:
                log_safe(self.app.state, "Python: No new junctions added (all already present).")
        except Exception as e:
            log_safe(self.app.state, f"ERROR in on_add_j_click: {e}")
            sys.stderr.flush()

    def on_remove_j_click(self, event=None):
        """Callback for 'Remove Selected' button in sidebar."""
        
        selected = self.L["mul_curr_j"].value
        log_safe(self.app.state, f"Removing indices: {selected}")
        if not selected: return
        new_list = []
        to_remove = set()
        for s in selected:
            parts = s.split('-')
            to_remove.add((int(parts[0]), int(parts[1])))
        for j in self.L["user_junctions"]:
            if j not in to_remove: new_list.append(j)
        self.L["user_junctions"] = new_list
        log_safe(self.app.state, f"Removed {len(to_remove)} junctions.")
        self.update_junction_ui()

    def on_mul_curr_j_change(self, attr, old, new):
        """Callback when the selection in the sidebar junction list changes."""
        if not new: return
        
        # 1. Zoom to the first selected junction (if enabled)
        if 0 in self.L["chk_autofocus_creator"].active:
            parts = new[0].split('-')
            target = (int(parts[0]), int(parts[1]))
            padding = (target[1] - target[0]) * 0.5
            self.L["shared_x_range"].start = target[0] - padding
            self.L["shared_x_range"].end = target[1] + padding

        # 2. Sync selection back to junction diagrams and creator highlighing
        ref_name = "unknown"
        for s_n, objs in self.L.get("sample_plots", {}).items():
            ds_js = objs[2]
            if ds_js.data and 'reference' in ds_js.data and len(ds_js.data['reference']) > 0:
                ref_name = ds_js.data['reference'][0]
                break
        
        junction_tuples = []
        creator_indices = []
        ds_creator = self.L["ds_creator_annotations"]
        creator_data = ds_creator.data
        
        for val in new:
            p = val.split('-')
            if len(p) != 2: continue
            js, je = int(p[0]), int(p[1])
            junction_tuples.append((ref_name, js, je))
            
            # Find index in ds_creator_annotations to highlight it there too
            if creator_data and 'left' in creator_data:
                for i in range(len(creator_data['left'])):
                    if creator_data['feature'][i] == 'Intron' and \
                       int(creator_data['left'][i]) == js and \
                       int(creator_data['right'][i]) == je:
                        creator_indices.append(i)
                        break

        # Highlight in creator plot if not already selected
        if list(ds_creator.selected.indices) != creator_indices:
            ds_creator.selected.indices = creator_indices
            
        # Update creator info div if we found a match
        if creator_indices:
            idx = creator_indices[0]
            text = f"<b>Selection:</b> {creator_data['feature'][idx]} {creator_data['feature_num'][idx]} | "
            text += f"Location: {int(creator_data['left'][idx]):,} - {int(creator_data['right'][idx]):,} | "
            text += f"Size: {int(creator_data['width'][idx]):,} bp"
            self.L["div_creator_info"].visible = True
            self.L["div_creator_info"].text = text

        # Highlight in coverage plots and potentially fetch reads
        if junction_tuples:
            self.app.on_junction_selection_change("creator", new=junction_tuples)

    def on_clear_j_click(self, event=None):
        """Clears all data in the Transcript Creator."""
        try:
            log_safe(self.app.state, "Python: Clearing Transcript Creator data...")
            

            # Clear state
            self.L["user_junctions"] = []
            self.L["num_t_start"].value = None
            self.L["num_t_end"].value = None
            
            if "div_creator_info" in self.L and self.L["div_creator_info"] is not None:
                self.L["div_creator_info"].visible = False
                self.L["div_creator_info"].text = ""
            
            self.update_junction_ui()

            # Clear selection in all annotation data sources
            for ds_name in ["ds_annotations", "ds_cds", "ds_transcripts", "ds_intron_markers"]:
                try:
                    self.L[ds_name].selected.indices = []
                except Exception:
                    pass

            log_safe(self.app.state, "Python: Transcript Creator cleared.")

        except Exception as e:
            log_safe(self.app.state, f"CRITICAL ERROR in clear: {e}")

    def on_create_t_click(self, event=None):
        """Finalizes the custom transcript and injects it into the global annotation state."""
        log_safe(self.app.state, "Button: Create and Inject Transcript clicked")
        
        if not self.L["user_junctions"]:
            log_safe(self.app.state, "Cannot create transcript: no junctions added.")
            return
        
        new_tid = self.L["txt_transcript_id"].value.strip()
        if not new_tid:
            log_safe(self.app.state, "Please provide a Transcript ID.")
            return

        junctions = sorted(self.L["user_junctions"])
        current_gene = self.L["sel_gene"].value
        current_gene_id = "User_Gene"
        strand = self.L["sel_strand"].value if "sel_strand" in self.L else "+"
        
        with self.app.state["lock"]:
            if self.app.state["processed_annotations"] is not None:
                df = self.app.state["processed_annotations"]
                mask = (df['gene_name'] == current_gene) if 'gene_name' in df.columns else (df['gene_id'] == current_gene)
                if any(mask):
                    sub = df[mask]
                    current_gene_id = sub['gene_id'].iloc[0] if 'gene_id' in sub.columns else current_gene
        
        new_annotations = []
        manual_start = self.L["num_t_start"].value
        manual_end = self.L["num_t_end"].value
        start_pos = manual_start if manual_start is not None else (junctions[0][0] - 50)
        end_pos = manual_end if manual_end is not None else (junctions[-1][1] + 50)
        
        if start_pos >= junctions[0][0]: start_pos = junctions[0][0] - 1
        if end_pos <= junctions[-1][1]: end_pos = junctions[-1][1] + 1

        new_annotations.append({
            "transcript_id": new_tid, "feature": "exon", "start": start_pos, "end": junctions[0][0],
            "strand": strand, "gene_id": current_gene_id, "gene_name": current_gene
        })
        for i in range(len(junctions) - 1):
            e_prev = junctions[i][1]; s_next = junctions[i+1][0]
            if s_next > e_prev:
                new_annotations.append({
                    "transcript_id": new_tid, "feature": "exon", "start": e_prev, "end": s_next,
                    "strand": strand, "gene_id": current_gene_id, "gene_name": current_gene
                })
        new_annotations.append({
            "transcript_id": new_tid, "feature": "exon", "start": junctions[-1][1], "end": end_pos,
            "strand": strand, "gene_id": current_gene_id, "gene_name": current_gene
        })
        
        with self.app.state["lock"]:
            if not self.app.state["analysis_data"]: self.app.state["analysis_data"] = {}
            if "annotations" not in self.app.state["analysis_data"]: self.app.state["analysis_data"]["annotations"] = []
            self.app.state["analysis_data"]["annotations"].extend(new_annotations)
            
        log_safe(self.app.state, f"Created transcript {new_tid} with {len(new_annotations)} exons. Refreshing...")
        process_analysis_data()
        self.app.update_all_samples()
        self.on_clear_j_click()

    def on_fetch_start_click(self, event=None):
        """Sets the transcript start to the left-most visible coordinate in the shared X range."""
        log_safe(self.app.state, "Button: Use View Start clicked")
        
        xr = self.L["shared_x_range"]
        if xr and xr.start != 0:
            self.L["num_t_start"].value = int(xr.start)
            log_safe(self.app.state, f"Start set from view: {int(xr.start)}")
            
    def on_fetch_end_click(self, event=None):
        """Sets the transcript end to the right-most visible coordinate in the shared X range."""
        log_safe(self.app.state, "Button: Use View End clicked")
        
        xr = self.L["shared_x_range"]
        if xr and xr.end != 1:
            self.L["num_t_end"].value = int(xr.end)
            log_safe(self.app.state, f"End set from view: {int(xr.end)}")
            
    def on_zoom_to_creator_click(self, event=None):
        """Zooms all plots to the bounds currently defined in the creator (start and end)."""
        log_safe(self.app.state, "Button: Zoom to Creator Bounds clicked")
        
        s = self.L["num_t_start"].value; e = self.L["num_t_end"].value
        if s is not None and e is not None and e > s:
            self.L["shared_x_range"].start = s; self.L["shared_x_range"].end = e
            log_safe(self.app.state, f"Zoomed view to creator bounds: {s} - {e}")

    def on_export_json_click(self, event=None):
        """Exports state for re-importing."""
        import json
        try:
            data = {
                "transcript_id": self.L["txt_transcript_id"].value,
                "strand": self.L["sel_strand"].value,
                "start": self.L["num_t_start"].value,
                "end": self.L["num_t_end"].value,
                "junctions": self.L["user_junctions"]
            }
            self.L["div_download"].text = f"json::CONTENT::{json.dumps(data)}"
        except Exception as e:
            log_safe(self.app.state, f"Export failed: {e}")

    def on_export_gtf_click(self, event=None):
        """Generates a GTF formatted string."""
        self._export_genomic_format("gtf")

    def on_export_gff3_click(self, event=None):
        """Generates a GFF3 formatted string."""
        self._export_genomic_format("gff3")

    def _export_genomic_format(self, fmt):
        try:
            tid = self.L["txt_transcript_id"].value.strip() or "New_Transcript"
            gene_name = self.L["sel_gene"].value or "User_Gene"
            strand = self.L["sel_strand"].value
            junctions = sorted(self.L["user_junctions"])
            
            if not junctions:
                log_safe(self.app.state, "Nothing to export: add junctions first.")
                return

            # Get reference/chromosome name
            ref = "unknown"
            for _, objs in self.L.get("sample_plots", {}).items():
                if objs[2].data.get('reference'):
                    ref = objs[2].data['reference'][0]; break

            m_start = self.L["num_t_start"].value
            m_end = self.L["num_t_end"].value
            start_pos = (m_start if m_start is not None else junctions[0][0] - 50)
            end_pos = (m_end if m_end is not None else junctions[-1][1] + 50)

            exons = []
            exons.append((start_pos, junctions[0][0]))
            for i in range(len(junctions) - 1):
                exons.append((junctions[i][1], junctions[i+1][0]))
            exons.append((junctions[-1][1], end_pos))

            lines = []
            if fmt == "gtf":
                for i, (s, e) in enumerate(exons):
                    attr = f'gene_id "{gene_name}"; transcript_id "{tid}"; exon_number "{i+1}";'
                    # GTF is 1-based inclusive. Internally we are 0-based inclusive, 0-based exclusive.
                    # 1-based [S, E] == 0-based [S-1, E) -> S = s+1, E = e
                    lines.append(f"{ref}\tlncSeeker\texon\t{int(s) + 1}\t{int(e)}\t.\t{strand}\t.\t{attr}")
            else: # GFF3
                lines.append("##gff-version 3")
                # GFF3 is also 1-based inclusive
                lines.append(f"{ref}\tlncSeeker\ttranscript\t{int(start_pos) + 1}\t{int(end_pos)}\t.\t{strand}\t.\tID={tid};Name={tid};Parent={gene_name}")
                for i, (s, e) in enumerate(exons):
                    lines.append(f"{ref}\tlncSeeker\texon\t{int(s) + 1}\t{int(e)}\t.\t{strand}\t.\tID={tid}.exon{i+1};Parent={tid}")

            self.L["div_download"].text = f"{fmt}::CONTENT::" + "\n".join(lines)
            log_safe(self.app.state, f"Generated {fmt.upper()} for {tid}")
        except Exception as e:
            log_safe(self.app.state, f"Genomic export failed: {e}")

    def on_import_json_change(self, attr, old, new):
        """Processes the uploaded JSON file and populates the creator state."""
        import json
        import base64
        if not new: return
        try:
            # Bokeh FileInput 'value' is base64 encoded string
            decoded = base64.b64decode(new).decode('utf-8')
            data = json.loads(decoded)
            
            # Basic validation
            if "junctions" not in data:
                log_safe(self.app.state, "Import failed: Invalid format (missing 'junctions').")
                return

            self.L["txt_transcript_id"].value = data.get("transcript_id", "Imported_Transcript")
            self.L["sel_strand"].value = data.get("strand", "+")
            self.L["num_t_start"].value = data.get("start")
            self.L["num_t_end"].value = data.get("end")
            
            # Convert junctions back to list of tuples if needed (JSON makes them lists)
            juncs = [tuple(j) for j in data.get("junctions", [])]
            self.L["user_junctions"] = juncs
            
            self.update_junction_ui()
            log_safe(self.app.state, f"Python: Successfully imported transcript: {data.get('transcript_id')}")
        except Exception as e:
            log_safe(self.app.state, f"Import failed: {e}")

    def on_ann_click(self, attr, old, new):
        """Callback when an element in the main annotation plot is clicked. Imports transcript data."""
        if not new: return
        curr_j = self.L.get("user_junctions")
        if curr_j: return # Don't overwrite if junctions already present? Or maybe we should?
            
        tid = None
        for ds_name in ["ds_annotations", "ds_cds", "ds_transcripts", "ds_intron_markers"]:
            ds = self.L[ds_name]
            # Bokeh 3.x uses indices list
            if ds.selected.indices == new:
                idx = new[0]
                if idx < len(ds.data['transcript']):
                    tid = ds.data['transcript'][idx]; break
            
        if not tid: return
            
        found_junctions = []
        t_start, t_end = None, None
        t_strand = "+"
        with self.app.state["lock"]:
            df = self.app.state.get("flat_annotations")
            if df is not None:
                t_df = df[df['transcript_id'] == tid].sort_values('start')
                starts = t_df['start'].tolist()
                ends = t_df['end'].tolist()
                if not t_df.empty: t_strand = t_df['strand'].iloc[0]
                if starts: t_start = int(min(starts)); t_end = int(max(ends))
                for i in range(len(starts) - 1):
                    found_junctions.append((int(ends[i]), int(starts[i+1])))

        if found_junctions:
            self.L["user_junctions"] = found_junctions
            if t_start is not None: self.L["num_t_start"].value = t_start; self.L["num_t_end"].value = t_end
            if "sel_strand" in self.L: self.L["sel_strand"].value = t_strand
            if "txt_transcript_id" in self.L: self.L["txt_transcript_id"].value = f"{tid}_new"
            self.update_junction_ui()
            log_safe(self.app.state, f"Imported {len(found_junctions)} junctions for transcript {tid}.")
        
        # Sync selection across annotation sources
        for ds_name in ["ds_annotations", "ds_cds", "ds_transcripts", "ds_intron_markers"]:
            ds = self.L[ds_name]
            if 'transcript' in ds.data:
                target_indices = [i for i, v in enumerate(ds.data['transcript']) if v == tid]
                if list(ds.selected.indices) != target_indices:
                    ds.selected.indices = target_indices

    def on_creator_ann_select(self, attr, old, new):
        """Handle selection of elements in the Transcript Creator preview plot."""
        if not new: return
        data = self.L["ds_creator_annotations"].data
        self.L["div_creator_info"].visible = True
        idx = new[0]
        text = f"<b>Selection:</b> {data['feature'][idx]} {data['feature_num'][idx]} | "
        text += f"Location: {int(data['left'][idx]):,} - {int(data['right'][idx]):,} | "
        text += f"Size: {int(data['width'][idx]):,} bp"
        self.L["div_creator_info"].text = text

        # Sync selection back to junction diagrams and list
        ref_name = "unknown"
        for s_n, objs in self.L.get("sample_plots", {}).items():
            ds_js = objs[2]
            if ds_js.data and 'reference' in ds_js.data and len(ds_js.data['reference']) > 0:
                ref_name = ds_js.data['reference'][0]
                break
        
        junction_tuples = []
        mul_select_vals = []
        for i in new:
            if i < len(data['feature']) and data['feature'][i] == 'Intron':
                js, je = int(data['left'][i]), int(data['right'][i])
                junction_tuples.append((ref_name, js, je))
                mul_select_vals.append(f"{js}-{je}")
        
        if junction_tuples:
            # Update the sidebar list
            self.L["mul_curr_j"].value = mul_select_vals
            # Highlight in coverage plots
            self.app.on_junction_selection_change("creator", new=junction_tuples)

    def update_preview_data(self, creator_state, creator_changed):
        """Updates the data sources for the creator preview plot."""
        import math
        import pandas as pd

        # We assume creator_changed has been determined by the app
        if creator_changed:
            self.L["last_creator_state"] = creator_state
            junctions = sorted(self.L.get("user_junctions", []))
            if junctions:
                manual_start = self.L["num_t_start"].value
                manual_end = self.L["num_t_end"].value
                strand = self.L["sel_strand"].value if "sel_strand" in self.L else "+"
                start_pos = manual_start if manual_start is not None else (junctions[0][0] - 50)
                end_pos = manual_end if manual_end is not None else (junctions[-1][1] + 50)
                if start_pos >= junctions[0][0]:
                    start_pos = junctions[0][0] - 1
                if end_pos <= junctions[-1][1]:
                    end_pos = junctions[-1][1] + 1
                
                creator_features = []
                creator_features.append({'start': start_pos, 'end': junctions[0][0], 'type': 'Exon'})
                for i in range(len(junctions)):
                    creator_features.append({'start': junctions[i][0], 'end': junctions[i][1], 'type': 'Intron'})
                    next_exon_start = junctions[i][1]
                    next_exon_end = junctions[i+1][0] if (i+1) < len(junctions) else end_pos
                    if next_exon_end > next_exon_start:
                        creator_features.append({'start': next_exon_start, 'end': next_exon_end, 'type': 'Exon'})
                
                tid = self.L["txt_transcript_id"].value.strip() or "New_Transcript"
                num_exons = sum(1 for f in creator_features if f['type'] == 'Exon')
                num_introns = sum(1 for f in creator_features if f['type'] == 'Intron')
                exon_counter = 1 if strand == '+' else num_exons
                intron_counter = 1 if strand == '+' else num_introns
                exon_step = 1 if strand == '+' else -1
                intron_step = 1 if strand == '+' else -1
                existing_js = set()
                for s_n in self.L["sel_samples"].value:
                    sample_d = self.app.selection_manager.get_sample_data(s_n)
                    if sample_d:
                        for js in sample_d.get("junction_spans", []):
                            existing_js.add((int(js['start']), int(js['end'])))
                
                lefts, rights, mids, widths, hs, ftrs, f_nums, colors, s_colors = [], [], [], [], [], [], [], [], []
                for f in creator_features:
                    lefts.append(f['start'])
                    rights.append(f['end'])
                    mids.append((f['start'] + f['end']) / 2)
                    widths.append(max(1, f['end'] - f['start']))
                    ftrs.append(f['type'])
                    if f['type'] == 'Exon':
                        f_nums.append(str(exon_counter))
                        exon_counter += exon_step
                        hs.append(0.25)
                        colors.append("#4CAF50")
                        s_colors.append("#4CAF50")
                    else:
                        f_nums.append(str(intron_counter))
                        intron_counter += intron_step
                        hs.append(0.04)
                        s_colors.append("gold")
                        j_key = (int(f['start']), int(f['end']))
                        if j_key in existing_js:
                            colors.append("gray")
                        else:
                            colors.append("salmon")
                
                self.L["ds_creator_annotations"].data = dict(left=lefts, right=rights, mid=mids, width=widths, h=hs, y=[0]*len(lefts), label=["Creator"]*len(lefts), transcript=[tid]*len(lefts), feature=ftrs, feature_num=f_nums, strand=[strand]*len(lefts), color=colors, sel_color=s_colors)
                self.L["ds_creator_transcripts"].data = dict(start=[min(start_pos, end_pos)], end=[max(start_pos, end_pos)], y=[0], strand=[strand], transcript=[tid])
                creator_markers = []
                m_angle = -math.pi/2 if strand == '+' else math.pi/2
                for i in range(len(junctions)):
                    j_start, j_end = junctions[i][0], junctions[i][1]
                    if j_end > j_start:
                        # Donor marker
                        creator_markers.append({'x': j_start, 'y': 0, 'angle': m_angle, 'marker': 'triangle_dot', 'size': 14})
                        # Acceptor marker
                        creator_markers.append({'x': j_end, 'y': 0, 'angle': m_angle, 'marker': 'inverted_triangle', 'size': 14})
                        
                        length = j_end - j_start
                        num_m = min(100, max(1, length // 2000))
                        step = length / (num_m + 1)
                        for m_idx in range(1, num_m + 1):
                            creator_markers.append({'x': j_start + m_idx * step, 'y': 0, 'angle': m_angle, 'marker': 'triangle', 'size': 8})
                self.L["ds_creator_markers"].data = dict(
                    x=[m['x'] for m in creator_markers], 
                    y=[m['y'] for m in creator_markers], 
                    angle=[m['angle'] for m in creator_markers], 
                    marker=[m['marker'] for m in creator_markers],
                    size=[m['size'] for m in creator_markers],
                    transcript=[tid]*len(creator_markers), 
                    strand=[strand]*len(creator_markers)
                )
            else:
                self.L["ds_creator_annotations"].data = dict(left=[], right=[], mid=[], width=[], h=[], y=[], label=[], transcript=[], feature_num=[], feature=[], strand=[], color=[], sel_color=[])
                self.L["ds_creator_transcripts"].data = dict(start=[], end=[], y=[], strand=[], transcript=[])
                self.L["ds_creator_markers"].data = dict(x=[], y=[], angle=[], marker=[], size=[], transcript=[], strand=[])

    def get_or_create_preview_plot(self, fixed_width, target_width, backend):

        """Ensures the creator preview plot exists and is configured correctly."""
        from bokeh.plotting import figure
        from bokeh.models import Range1d, NumeralTickFormatter, PanTool, WheelZoomTool, CustomJS, TapTool, HoverTool
        from bokeh.events import Reset
        from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot
        import math

        if self.L.get("p_creator") is None:
            p_creator = figure(height=120, sizing_mode="fixed" if fixed_width else "stretch_width", 
                               width=target_width if fixed_width else None, title="Transcript Editor Preview", 
                               y_range=Range1d(start=-0.5, end=0.5), 
                               x_range=self.L["shared_x_range"], tools="reset,save,tap", output_backend=backend)
            self.L["p_creator"] = p_creator
            p_creator.xaxis.formatter = NumeralTickFormatter(format="0")
            p_creator.xaxis.major_label_orientation = math.pi/4
            p_creator.xgrid.grid_line_color = None
            p_creator.ygrid.grid_line_color = None
            p_creator.yaxis.visible = False
            
            pan_creator = PanTool(description="Pan")
            w_zoom_creator = WheelZoomTool(dimensions="width")
            p_creator.add_tools(pan_creator, w_zoom_creator)
            p_creator.toolbar.active_drag = pan_creator
            p_creator.toolbar.active_scroll = w_zoom_creator
            
            p_creator.segment(x0='start', y0='y', x1='end', y1='y', source=self.L["ds_creator_transcripts"], color="gray")
            self.L["r_creator_markers"] = p_creator.scatter(x='x', y='y', size='size', marker='marker', angle='angle', source=self.L["ds_creator_markers"], color="gray", alpha=0.8)
            self.L["r_creator_ann"] = p_creator.rect(x='mid', y='y', width='width', height='h', source=self.L["ds_creator_annotations"], color="color", alpha=0.7, selection_color="sel_color", nonselection_color="color", selection_alpha=0.9, nonselection_alpha=0.7)
            
            tap_creator = TapTool(renderers=[self.L["r_creator_ann"]])
            p_creator.add_tools(tap_creator)
            
            hover_c = HoverTool(renderers=[self.L["r_creator_ann"]], tooltips=[
                ("Transcript", "@transcript"), ("Type", "@feature"), ("Number", "@feature_num"), 
                ("Strand", "@strand"), ("Location", "@left{0,0} - @right{0,0}"), ("Size", "@width{0,0} bp")
            ])
            p_creator.add_tools(hover_c)
            add_crosshair_to_plot(p_creator, self.L, self.app.state, plot_type="creator")
            
            self.L["ds_creator_annotations"].selected.on_change('indices', self.on_creator_ann_select)
            p_creator.on_event(Reset, self.app.on_reset_zoom_click)
            p_creator.js_on_event(Reset, CustomJS(args=dict(xr=p_creator.x_range, dr=self.L["ds_core_range"]), code="""if (dr.data['start'].length > 0) { xr.start = dr.data['start'][0]; xr.end = dr.data['end'][0]; }"""))
        else:
            p_creator = self.L["p_creator"]
            if p_creator.x_range != self.L["shared_x_range"]:
                p_creator.x_range = self.L["shared_x_range"]
            p_creator.output_backend = backend
            if fixed_width:
                p_creator.sizing_mode = "fixed"
                p_creator.width = target_width
            else:
                p_creator.sizing_mode = "stretch_width"
            
            if "r_creator_ann" in self.L:
                for tool in p_creator.tools:
                    if isinstance(tool, TapTool): 
                        tool.renderers = [self.L["r_creator_ann"]]
        
        if "chk_xaxis_creator" in self.L: 
            p_creator.xaxis.visible = (0 in self.L["chk_xaxis_creator"].active)
        if "chk_markers_creator" in self.L and "r_creator_markers" in self.L: 
            self.L["r_creator_markers"].visible = (0 in self.L["chk_markers_creator"].active)
            
        return p_creator

