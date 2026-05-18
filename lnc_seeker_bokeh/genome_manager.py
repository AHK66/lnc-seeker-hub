# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import math
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import (
    Range1d, NumeralTickFormatter, PanTool, WheelZoomTool, 
    CustomJS, LabelSet, HoverTool, CustomAction
)
from bokeh.events import Reset, MouseWheel
from lnc_seeker_bokeh.state import log_safe, strip_id
from lnc_seeker_bokeh.data_utils import process_analysis_data
from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot, get_transcript_color
from lnc_seeker_bokeh.genome_plot import create_genome_annotations_plot

class GenomeManager:
    def __init__(self, app):
        """
        Coordinates the 'Genome Annotations' plot.
        :param app: The VisualizerApp instance this manager belongs to.
        """
        self.app = app
        self.L = app.L

    def on_ann_filter_change(self, attr, old, new):
        """Handles changes in annotation visibility filtering."""
        active = new
        self.app.state["config"]["genome_annotations"]["show_full_range"] = (0 in active)
        process_analysis_data(self.app.state)
        with self.app.state["lock"]:
            if self.app.state["processed_annotations"] is not None:
                df = self.app.state["processed_annotations"]
                if not df.empty:
                    min_x, max_x = df['start'].min(), df['end'].max()
                    self.L["shared_x_range"].start = min_x
                    self.L["shared_x_range"].end = max_x
            self.L["reset_y_range_ann"] = True
        self.app.update_all_samples()

    def on_reset_ann_click(self):
        """Resets x-range and annotation y-range."""
        with self.app.state["lock"]:
            data = self.app.state.get("analysis_data", {})
            min_x = data.get("min_x")
            max_x = data.get("max_x")
            
            if min_x is not None and max_x is not None and max_x > min_x:
                if abs(self.L["shared_x_range"].start - min_x) > 1.0 or abs(self.L["shared_x_range"].end - max_x) > 1.0:
                    log_safe(self.app.state, f"Resetting view to core area: {min_x}..{max_x}")
                    self.L["shared_x_range"].start = min_x
                    self.L["shared_x_range"].end = max_x
                
            # Reset annotation Y range
            if self.L.get("p_ann") and self.L.get("last_max_y_ann") is not None:
                self.L["p_ann"].y_range.start = -1
                self.L["p_ann"].y_range.end = self.L["last_max_y_ann"] + 1

    def update_genome_data(self, ann_changed, data_changed=False):
        with self.app.state["lock"]:
            df_exons = self.app.state["processed_annotations"]
            df_cds = self.app.state.get("processed_cds")
            df_t = self.app.state["processed_transcripts"]
            df_gl = self.app.state.get("processed_gene_labels")
            df_markers = self.app.state.get("processed_markers")
        
        if ann_changed:
            if df_exons is not None:
                # Filter out 'transcript' features entirely from the box/rect rendering.
                # Exons, UTRs, and ncRNA segments are drawn as blocks; the transcript 
                # backbone is drawn separately as a thin line.
                df_boxes = df_exons[df_exons['feature'] != 'transcript'].copy()

                if 'gene_name' in df_boxes.columns:
                    labels = df_boxes['gene_name'].fillna(df_boxes['gene_id']).fillna("unknown")
                else:
                    labels = df_boxes.get('gene_id', pd.Series(["unknown"] * len(df_boxes)))
                
                ex_nums = df_boxes['exon_number'].fillna("-") if 'exon_number' in df_boxes.columns else pd.Series(["-"] * len(df_boxes))
                colors = [get_transcript_color(self.app.state, tid) for tid in df_boxes['transcript_id']]
                
                # Configure height based on feature type
                f_alpha = self.app.state.get("config", {}).get("genome_annotations", {}).get("feature_alpha", 0.6)
                
                alphas = [f_alpha] * len(df_boxes)
                heights = [0.35] * len(df_boxes) # Standard height for exons/features

                self.L["ds_annotations"].data = dict(
                    left=df_boxes['start'].tolist(), right=df_boxes['end'].tolist(),
                    mid=((df_boxes['start'] + df_boxes['end']) / 2).tolist(),
                    width=(df_boxes['end'] - df_boxes['start']).tolist(),
                    y=df_boxes['y'].tolist(), label=labels.tolist(),
                    transcript=df_boxes['transcript_id'].tolist(),
                    exon_num=ex_nums.tolist(),
                    cds_num=pd.Series(["-"] * len(df_boxes)).tolist(),
                    feature=df_boxes['feature'].tolist(),
                    strand=df_boxes['strand'].tolist(),
                    color=colors,
                    alpha=alphas,
                    height=heights
                )
            else:
                self.L["ds_annotations"].data = dict(
                    left=[], right=[], mid=[], width=[], y=[], label=[], 
                    transcript=[], exon_num=[], cds_num=[], feature=[], 
                    strand=[], color=[], alpha=[], height=[]
                )

            if df_cds is not None:
                if 'gene_name' in df_cds.columns:
                    cds_labels = df_cds['gene_name'].fillna(df_cds['gene_id']).fillna("unknown")
                else:
                    cds_labels = df_cds.get('gene_id', pd.Series(["unknown"] * len(df_cds)))
                
                cds_ex_nums = df_cds['exon_number'].fillna("-") if 'exon_number' in df_cds.columns else pd.Series(["-"] * len(df_cds))
                cds_ranks = df_cds['cds_rank'] if 'cds_rank' in df_cds.columns else pd.Series(["-"] * len(df_cds))
                cds_colors = [get_transcript_color(self.app.state, tid) for tid in df_cds['transcript_id']] if 'transcript_id' in df_cds.columns else []
                
                self.L["ds_cds"].data = dict(
                    left=df_cds['start'].tolist() if 'start' in df_cds.columns else [],
                    right=df_cds['end'].tolist() if 'end' in df_cds.columns else [],
                    mid=((df_cds['start'] + df_cds['end']) / 2).tolist() if 'start' in df_cds.columns else [],
                    width=(df_cds['end'] - df_cds['start']).tolist() if 'start' in df_cds.columns else [],
                    y=df_cds['y'].tolist() if 'y' in df_cds.columns else [],
                    label=cds_labels.tolist(),
                    transcript=df_cds['transcript_id'].tolist() if 'transcript_id' in df_cds.columns else [],
                    exon_num=cds_ex_nums.tolist(),
                    cds_num=cds_ranks.tolist(),
                    feature=df_cds['feature'].tolist() if 'feature' in df_cds.columns else [],
                    strand=df_cds['strand'].tolist() if 'strand' in df_cds.columns else [],
                    color=cds_colors
                )
            else:
                self.L["ds_cds"].data = dict(left=[], right=[], mid=[], width=[], y=[], label=[], transcript=[], exon_num=[], cds_num=[], feature=[], strand=[], color=[])
            
            if df_t is not None:
                self.L["ds_transcripts"].data = dict(
                    start=df_t['start'].tolist() if 'start' in df_t.columns else [],
                    end=df_t['end'].tolist() if 'end' in df_t.columns else [],
                    y=df_t['y'].tolist() if 'y' in df_t.columns else [],
                    strand=df_t['strand'].tolist() if 'strand' in df_t.columns else [],
                    transcript=df_t['transcript_id'].tolist() if 'transcript_id' in df_t.columns else []
                )
            else:
                self.L["ds_transcripts"].data = dict(start=[], end=[], y=[], strand=[], transcript=[])

            if df_gl is not None:
                self.L["ds_gene_labels"].data = dict(
                    x=df_gl['x'].tolist() if 'x' in df_gl.columns else [],
                    y=df_gl['y'].tolist() if 'y' in df_gl.columns else [],
                    text=df_gl['text'].tolist() if 'text' in df_gl.columns else []
                )
            else:
                self.L["ds_gene_labels"].data = dict(x=[], y=[], text=[])

            if df_markers is not None:
                self.L["ds_intron_markers"].data = dict(
                    x=df_markers['x'].tolist() if 'x' in df_markers.columns else [],
                    y=df_markers['y'].tolist() if 'y' in df_markers.columns else [],
                    angle=df_markers['angle'].tolist() if 'angle' in df_markers.columns else [],
                    transcript=df_markers['transcript_id'].tolist() if 'transcript_id' in df_markers.columns else [],
                    strand=df_markers['strand'].tolist() if 'strand' in df_markers.columns else []
                )
            else:
                self.L["ds_intron_markers"].data = dict(x=[], y=[], angle=[], transcript=[], strand=[])

            # Update Intron segments for hover
            df_introns = self.app.state.get("processed_introns")
            if df_introns is not None:
                self.L["ds_introns"].data = dict(
                    left=df_introns['start'].tolist(),
                    right=df_introns['end'].tolist(),
                    mid=((df_introns['start'] + df_introns['end']) / 2).tolist(),
                    width=(df_introns['end'] - df_introns['start']).tolist(),
                    y=df_introns['y'].tolist(),
                    transcript=df_introns['transcript_id'].tolist(),
                    gene_id=df_introns['gene_id'].tolist(),
                    strand=df_introns['strand'].tolist(),
                    length=df_introns['length'].tolist(),
                    reference=df_introns['reference'].tolist()
                )
            else:
                self.L["ds_introns"].data = dict(left=[], right=[], mid=[], width=[], y=[], transcript=[], gene_id=[], strand=[], length=[], reference=[])

        if ann_changed or data_changed:
            self._update_unsupported_introns_data()

    def _update_unsupported_introns_data(self):
        """Identifies introns that lack matching junction evidence in current samples."""
        if "chk_unsupported_introns" not in self.L:
            return
            
        is_active = (0 in self.L["chk_unsupported_introns"].active)
        if not is_active:
            self.L["ds_unsupported_introns"].data = dict(x=[], y=[])
            return

        # 1. Gather all active junctions from analysis_data
        selected_samples = self.app.selection_manager.get_selected_samples()
        min_reads = self.L["sld_min_reads"].value
        use_clean = (0 in self.L["chk_filter_flanks"].active) if "chk_filter_flanks" in self.L else False
        
        analysis_data = self.app.state.get("analysis_data", {})
        samples_data = analysis_data.get("samples", {})
        
        active_junctions = set()
        for name in selected_samples:
            raw_name = strip_id(name)
            s_data = samples_data.get(raw_name)
            if not s_data:
                continue
            
            for js in s_data.get("junction_spans", []):
                reads = js.get('reads_clean', js['reads']) if use_clean else js['reads']
                if reads >= min_reads:
                    active_junctions.add((int(js['start']), int(js['end'])))
        
        # 2. Compare with introns
        df_introns = self.app.state.get("processed_introns")
        if df_introns is not None and not df_introns.empty:
            # An intron is unsupported if its (start, end) is NOT in active_junctions
            unsupported_mask = df_introns.apply(
                lambda row: (int(row['start']), int(row['end'])) not in active_junctions, 
                axis=1
            )
            df_unsupported = df_introns[unsupported_mask]
            
            self.L["ds_unsupported_introns"].data = dict(
                x=((df_unsupported['start'] + df_unsupported['end']) / 2).tolist(),
                y=df_unsupported['y'].tolist()
            )
        else:
            self.L["ds_unsupported_introns"].data = dict(x=[], y=[])

    def get_or_create_genome_plot(self, ann_target_h, fixed_width, target_width, backend):
        """Ensures the annotation plot exists and is configured correctly."""
        if self.L.get("p_ann") is None:
            p_ann = create_genome_annotations_plot(
                self.L, self.app.state, self.L["shared_x_range"], self.on_reset_ann_click, 
                fixed_width, target_width, ann_target_h, backend
            )
            self.L["p_ann"] = p_ann
            
            # Selection callback for introns
            self.L["ds_introns"].selected.on_change('indices', self.app.selection_manager.on_intron_selection_change)
        else:
            p_ann = self.L["p_ann"]

            if p_ann.x_range != self.L["shared_x_range"]:
                p_ann.x_range = self.L["shared_x_range"]
            p_ann.height = ann_target_h
            p_ann.output_backend = backend

            if fixed_width:
                p_ann.sizing_mode = "fixed"; p_ann.width = target_width
            else:
                p_ann.sizing_mode = "stretch_width"
        
        # Sync visibility and ranges
        try:
            if 'y' in self.L["ds_annotations"].data and len(self.L["ds_annotations"].data['y']) > 0:
                new_max_y = float(np.max(self.L["ds_annotations"].data['y']))
                if self.L.get("last_max_y_ann") is None or self.L.get("reset_y_range_ann") or abs(self.L["last_max_y_ann"] - new_max_y) > 0.01:
                    log_safe(self.app.state, f"Updating Genome Annotations Y-range to {new_max_y}")
                    p_ann.y_range.start = -1
                    p_ann.y_range.end = new_max_y + 1
                    self.L["last_max_y_ann"] = new_max_y
                    self.L["reset_y_range_ann"] = False
        except Exception:
            pass
            
        p_ann.yaxis.visible = False
        if "r_ann_markers" in self.L:
            self.L["r_ann_markers"].visible = (0 in self.L["chk_markers"].active) if "chk_markers" in self.L else False
        if "r_ann_cds" in self.L:
            self.L["r_ann_cds"].visible = (0 in self.L["chk_cds"].active) if "chk_cds" in self.L else True
        if "chk_xaxis_ann" in self.L:
            p_ann.xaxis.visible = (0 in self.L["chk_xaxis_ann"].active)
        
        return p_ann
        
        # Sync visibility and ranges
        try:
            if 'y' in self.L["ds_annotations"].data and len(self.L["ds_annotations"].data['y']) > 0:
                new_max_y = float(np.max(self.L["ds_annotations"].data['y']))
                if self.L.get("last_max_y_ann") is None or self.L.get("reset_y_range_ann") or abs(self.L["last_max_y_ann"] - new_max_y) > 0.01:
                    log_safe(self.app.state, f"Updating Genome Annotations Y-range to {new_max_y}")
                    p_ann.y_range.start = -1
                    p_ann.y_range.end = new_max_y + 1
                    self.L["last_max_y_ann"] = new_max_y
                    self.L["reset_y_range_ann"] = False
        except Exception:
            pass
            
        p_ann.yaxis.visible = False
        if "r_ann_markers" in self.L:
            self.L["r_ann_markers"].visible = (0 in self.L["chk_markers"].active) if "chk_markers" in self.L else False
        if "r_ann_cds" in self.L:
            self.L["r_ann_cds"].visible = (0 in self.L["chk_cds"].active) if "chk_cds" in self.L else True
        if "chk_xaxis_ann" in self.L:
            p_ann.xaxis.visible = (0 in self.L["chk_xaxis_ann"].active)
        
        return p_ann
