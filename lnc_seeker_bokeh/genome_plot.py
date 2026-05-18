# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import math
import numpy as np
from bokeh.plotting import figure
from bokeh.models import (
    Range1d, NumeralTickFormatter, PanTool, WheelZoomTool, 
    CustomJS, LabelSet, HoverTool,
    CustomAction
)
from bokeh.events import Reset, MouseWheel
from lnc_seeker_bokeh.state import log_safe
from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot

def create_genome_annotations_plot(L, state, shared_x_range, on_reset_ann_click, fixed_width, target_width, ann_target_h, backend):
    max_y = 1
    try:
        if 'y' in L["ds_annotations"].data and len(L["ds_annotations"].data['y']) > 0:
            max_y = float(np.max(L["ds_annotations"].data['y']))
    except:
        pass
        
    p_ann = figure(height=ann_target_h, sizing_mode="fixed" if fixed_width else "stretch_width", 
                   width=target_width if fixed_width else None,
                   title="Genome Annotations", 
                   y_range=Range1d(start=-1, end=max_y + 1), 
                   x_range=shared_x_range,
                   tools="reset,save,tap",
                   output_backend=backend)
    
    p_ann.xaxis.formatter = NumeralTickFormatter(format="0")
    p_ann.xaxis.major_label_orientation = math.pi/4
    p_ann.xgrid.grid_line_color = None
    p_ann.ygrid.grid_line_color = None
    
    pan = PanTool(description="Pan")
    w_zoom_dummy_ann = WheelZoomTool(speed=0)

    # Height Adjustment Tools: Swapped logic (Up=Decrease, Down=Increase)
    dec_h = CustomAction(description="Decrease Diagram Height", icon="chevron_up",
                         callback=CustomJS(args=dict(spn=L["spn_ann_height"]), code="""
        if (spn.value - 50 >= spn.low) spn.value -= 50;
    """))
    inc_h = CustomAction(description="Increase Diagram Height", icon="chevron_down",
                         callback=CustomJS(args=dict(spn=L["spn_ann_height"]), code="""
        if (spn.value + 50 <= spn.high) spn.value += 50;
    """))

    p_ann.add_tools(pan, w_zoom_dummy_ann, dec_h, inc_h)
    p_ann.on_event(Reset, on_reset_ann_click)
    
    p_ann.js_on_event(Reset, CustomJS(args=dict(xr=p_ann.x_range, dr=L["ds_core_range"]), code="""
        if (dr.data['start'].length > 0) {
            xr.start = dr.data['start'][0];
            xr.end = dr.data['end'][0];
        }
    """))
    
    zoom_pan_ann_js = CustomJS(args=dict(xr=p_ann.x_range, yr=p_ann.y_range), code="""
        const delta = cb_obj.delta;
        const shift = !!(cb_obj.modifiers && cb_obj.modifiers.shift);
        if (shift) {
            const factor = delta > 0 ? 0.90 : 1.10;
            const y = cb_obj.y;
            if (y != null) {
                const s = yr.start; const e = yr.end;
                yr.start = y - (y - s) * factor;
                yr.end = y + (e - y) * factor;
            }
        } else {
            const factor = delta > 0 ? 0.90 : 1.10;
            const x = cb_obj.x;
            if (x != null) {
                const s = xr.start; const e = xr.end;
                xr.start = x - (x - s) * factor;
                xr.end = x + (e - x) * factor;
            }
        }
    """)
    p_ann.js_on_event(MouseWheel, zoom_pan_ann_js)
    p_ann.toolbar.active_scroll = w_zoom_dummy_ann
    p_ann.toolbar.active_drag = pan
    
    # Selection styling for introns from config
    sel_cfg = state.get("config", {}).get("genome_annotations", {}).get("selection_style", {})
    sel_color = sel_cfg.get("color", "red")
    sel_alpha = sel_cfg.get("alpha", 0.6)
    sel_lcolor = sel_cfg.get("line_color", "firebrick")
    sel_lwidth = sel_cfg.get("line_width", 2)
    sel_height = sel_cfg.get("height", 0.3)

    # Transcript spine (central line)
    L["r_ann_transcripts"] = p_ann.segment(x0='start', y0='y', x1='end', y1='y', 
                                              source=L["ds_transcripts"], color="black", line_width=1)
    
    L["r_ann_markers"] = p_ann.scatter(x='x', y='y', size=5, marker='triangle', angle='angle', 
                                         source=L["ds_intron_markers"], color="black", alpha=0.8)
    
    # Intron Segment for hover and selection (Invisible by default, Red when selected)
    L["r_ann_intron_seg"] = p_ann.rect(x='mid', y='y', width='width', height=sel_height, 
                                            source=L["ds_introns"], color="black", alpha=0.0,
                                            selection_fill_color=sel_color, selection_fill_alpha=sel_alpha,
                                            selection_line_color=sel_lcolor, selection_line_width=sel_lwidth,
                                            selection_line_alpha=min(1.0, sel_alpha + 0.2),
                                            nonselection_fill_alpha=0.0, nonselection_line_alpha=0.0)

    # Exons / Transcript-boxes (if detailed exons missing)
    L["r_ann_exons"] = p_ann.rect(x='mid', y='y', width='width', height='height', 
                                      source=L["ds_annotations"], color="color", 
                                      fill_alpha='alpha', line_color="black", line_width=0.5, line_alpha=0.4)
    
    # CDS (Coding Sequence) - Thickest
    L["r_ann_cds"] = p_ann.rect(x='mid', y='y', width='width', height=0.5, 
                                    source=L["ds_cds"], color="color", 
                                    fill_alpha=0.8, line_color="black", line_width=1.0)
    
    # Unsupported Intron Markers (Red X)
    L["r_ann_unsupported"] = p_ann.scatter(x='x', y='y', size=12, marker="x", 
                                                source=L["ds_unsupported_introns"], 
                                                color="red", line_width=3)

    tooltips_ann = [
        ("Gene", "@label"), ("Transcript", "@transcript"), ("Exon (Overall)", "@exon_num"),
        ("Exon (Coding/CDS)", "@cds_num"), ("Feature", "@feature"), ("Strand", "@strand"),
        ("Location", "@left{0,0} - @right{0,0}")
    ]
    p_ann.add_tools(HoverTool(renderers=[L["r_ann_exons"], L["r_ann_cds"]], tooltips=tooltips_ann))
    
    # Separate hover tools for markers to avoid column mismatch errors
    p_ann.add_tools(HoverTool(
        renderers=[L["r_ann_markers"]],
        tooltips=[("Transcript", "@transcript"), ("Strand", "@strand"), ("Position", "@x{0,0}")]
    ))

    p_ann.add_tools(HoverTool(
        renderers=[L["r_ann_intron_seg"]],
        tooltips=[("Feature", "Intron"), ("Gene", "@gene_id"), ("Transcript", "@transcript"), ("Length", "@length bp"), ("Strand", "@strand"), ("Location", "@left{0,0} - @right{0,0}")]
    ))

    p_ann.add_layout(LabelSet(x='x', y='y', text='text', source=L["ds_gene_labels"],
                              text_font_size="9pt", text_font_style="bold",
                              x_offset=0, text_align="center", text_baseline="bottom"))
    
    # Extension markers (dictionary-based selection)
    if "ds_extension" in L:
        # We use constants for y0 and y1 to ensure they span the full vertical range 
        # of the linear plot (starting from bottom) without relying on log-safe data columns
        r_ext_ann = p_ann.segment(x0='x0', y0=-1000, x1='x1', y1=1000000, source=L["ds_extension"],
                                  line_color='color', line_width='width', line_dash='dash', alpha='alpha',
                                  name="ext_segment")
        p_ann.add_tools(HoverTool(renderers=[r_ext_ann], tooltips=[
            ("Boundary", "@desc"),
            ("Extraction Region", "@region"),
            ("Position", "@x0{0,0}")
        ]))

    add_crosshair_to_plot(p_ann, L, state, plot_type="genome")
    return p_ann
