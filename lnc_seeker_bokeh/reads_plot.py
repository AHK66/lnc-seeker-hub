# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import math
from bokeh.plotting import figure
from bokeh.models import (
    Range1d, PanTool, WheelZoomTool, CustomJS, 
    NumeralTickFormatter, LabelSet, CDSView, CustomJSFilter, GroupFilter,
    CustomAction
)
from bokeh.events import Reset, MouseWheel
from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot

def create_reads_plot(L, state, reads_target_h, fixed_width, target_width, backend, on_reset_callback):
    """Initializes and configures the 'Full Read Layout' figure."""
    cfg_ru = state["config"]["full_read_layout"]
    lw_ru = cfg_ru.get("show_reads_line_width", 2)
    gap_ru = cfg_ru.get("vertical_gap", 4)
    squeeze_ru = cfg_ru.get("reads_vertical_squeeze", 0.2)
    
    # Calculate initial visible rows based on squeeze factor
    initial_visible_rows = (reads_target_h / (lw_ru + gap_ru)) / squeeze_ru
    
    p_reads = figure(height=reads_target_h, sizing_mode="fixed" if fixed_width else "stretch_width", 
                     width=target_width if fixed_width else None, title="Click a junction label to see full read layout", 
                     x_range=L["shared_x_range"], y_range=Range1d(start=0, end=initial_visible_rows), 
                     output_backend=backend, tools="reset,save,hover")
    
    pan_ru = PanTool(dimensions="both", description="Pan (x & y)")
    w_z_d_ru = WheelZoomTool(speed=0)
    
    # Vertical zoom (squeezing) logic with limits from config
    min_sq = cfg_ru.get("min_vertical_squeeze", 0.05)
    max_sq = cfg_ru.get("max_vertical_squeeze", 1.2)
    
    p_reads.js_on_event(MouseWheel, CustomJS(
        args=dict(xr=p_reads.x_range, yr=p_reads.y_range, busy=L["ds_reads_busy"], 
                  p=p_reads, lw=lw_ru, gap=gap_ru, min_sq=min_sq, max_sq=max_sq), 
        code="""
        if (busy.data['is_busy'][0] === 1) return;
        const delta = cb_obj.delta; 
        const shift = !!(cb_obj.modifiers && cb_obj.modifiers.shift); 
        if (shift) { 
            const unit = lw + gap;
            const h = p.height;
            const min_span = (h / unit) / max_sq;
            const max_span = (h / unit) / min_sq;
            
            const factor = delta > 0 ? 0.90 : 1.10; 
            const y = cb_obj.y; 
            if (y != null) { 
                const s = yr.start; const e = yr.end; 
                let new_s = y - (y - s) * factor; 
                let new_e = y + (e - y) * factor;
                const new_span = new_e - new_s;
                
                if (new_span < min_span || new_span > max_span) {
                    const clamped_span = Math.max(min_span, Math.min(max_span, new_span));
                    const ratio = clamped_span / (e - s);
                    new_s = y - (y - s) * ratio;
                    new_e = y + (e - y) * ratio;
                }
                yr.start = new_s; yr.end = new_e;
            } 
        } else { 
            const factor = delta > 0 ? 0.90 : 1.10; 
            const x = cb_obj.x; 
            if (x != null) { 
                const s = xr.start; const e = xr.end; 
                xr.start = x - (x - s) * factor; xr.end = x + (e - x) * factor; 
            } 
        }
    """))

    # Height Adjustment Tools: Swapped logic (Up=Decrease, Down=Increase)
    dec_h_ru = CustomAction(description="Decrease Diagram Height", icon="chevron_up",
                            callback=CustomJS(args=dict(spn=L["spn_reads_height"]), code="""
        if (spn.value - 50 >= spn.low) spn.value -= 50;
    """))
    inc_h_ru = CustomAction(description="Increase Diagram Height", icon="chevron_down",
                            callback=CustomJS(args=dict(spn=L["spn_reads_height"]), code="""
        if (spn.value + 50 <= spn.high) spn.value += 50;
    """))

    p_reads.add_tools(pan_ru, w_z_d_ru, dec_h_ru, inc_h_ru)
    p_reads.toolbar.active_scroll = w_z_d_ru
    p_reads.toolbar.active_drag = pan_ru
    p_reads.on_event(Reset, on_reset_callback)
    p_reads.js_on_event(Reset, CustomJS(args=dict(xr=p_reads.x_range, dr=L["ds_core_range"]), code="""
        if (dr.data['start'].length > 0) { xr.start = dr.data['start'][0]; xr.end = dr.data['end'][0]; }
    """))
    
    r_styles_ru = state["config"]["full_read_layout"].get("layout_styles", {})
    L["r_reads_exons"] = p_reads.segment(x0='x0', y0='y', x1='x1', y1='y', source=L["ds_reads_exons"], 
                                             color='color', line_width='thickness', alpha=r_styles_ru.get("exon_single", {}).get("alpha", 0.6), 
                                             legend_label="Aligned Segments")
    
    it_dash_ru = r_styles_ru.get("intron_target", {}).get("line_dash", "solid")
    if isinstance(it_dash_ru, list): it_dash_ru = tuple(it_dash_ru)
    
    L["r_reads_introns"] = p_reads.segment(x0='x0', y0='y', x1='x1', y1='y', source=L["ds_reads_introns"], 
                                               color='color', line_width='line_w', line_dash="solid", alpha=0.6, 
                                               legend_label="Within-Read Gaps")
    L["r_reads_introns"].glyph.line_dash = it_dash_ru

    # Deletion markers (red symbol)
    del_s = r_styles_ru.get("deletion", {})
    del_view = CDSView(filter=GroupFilter(column_name="type", group="deletion"))
    
    # We use scatter for the red symbol
    marker_type = del_s.get("marker", "x")
    L["r_reads_deletion_markers"] = p_reads.scatter(x='mid_x', y='y', source=L["ds_reads_introns"], view=del_view,
                                                   marker=marker_type, size=14,
                                                   line_color=del_s.get("color", "red"), 
                                                   fill_color=del_s.get("color", "red"),
                                                   line_width=3, alpha=1.0, 
                                                   legend_label="Deletions")
    
    p_reads.add_layout(LabelSet(x='mid_x', y='y', text='label', source=L["ds_reads_introns"], 
                                text_font_size=r_styles_ru.get("label_intron", {}).get("font_size", "7pt"), 
                                text_color=r_styles_ru.get("label_intron", {}).get("color", "black"), 
                                text_align="center", text_baseline="bottom"))
    
    # Mismatch / Insertion labels on exons (Colored box with white text)
    L["r_reads_exons_labels"] = LabelSet(x='mid_x', y='y', text='label', source=L["ds_reads_exons"],
                                         text_font_size=r_styles_ru.get("label_intron", {}).get("font_size", "7pt"),
                                         text_color="white", text_align="center", text_baseline="middle",
                                         background_fill_color='label_bg_color', background_fill_alpha=1.0)
    p_reads.add_layout(L["r_reads_exons_labels"])
    
    mb_s_ru = r_styles_ru.get("mate_bridge", {})
    mb_dash_ru = mb_s_ru.get("line_dash", "solid")
    if isinstance(mb_dash_ru, list): mb_dash_ru = tuple(mb_dash_ru)
    
    L["r_reads_bridges"] = p_reads.segment(x0='x0', y0='y', x1='x1', y1='y', source=L["ds_reads_bridges"], 
                                               color=mb_s_ru.get("color", "red"), line_width=mb_s_ru.get("line_width", 1), 
                                               alpha=mb_s_ru.get("alpha", 0.4), line_dash="solid", 
                                               legend_label="Mate Bridges")
    L["r_reads_bridges"].glyph.line_dash = mb_dash_ru
    
    L["r_reads_bridges_labels"] = LabelSet(x='mid_x', y='y', text='label', source=L["ds_reads_bridges"], 
                                                text_font_size=r_styles_ru.get("label_bridge", {}).get("font_size", "7pt"), 
                                                text_color=r_styles_ru.get("label_bridge", {}).get("color", "red"), 
                                                text_align="center", text_baseline="bottom")
    p_reads.add_layout(L["r_reads_bridges_labels"])
    
    sl_s_ru = r_styles_ru.get("label_sample", {})
    L["r_reads_labels"] = LabelSet(x='x', y='y', text='text', source=L["ds_reads_labels"], 
                                        text_font_size=sl_s_ru.get("font_size", "10pt"), text_font_style="bold", 
                                        text_color=sl_s_ru.get("color", "white"), text_baseline="bottom",
                                        background_fill_color=sl_s_ru.get("background_color", "black"), 
                                        background_fill_alpha=sl_s_ru.get("background_alpha", 1.0))
    p_reads.add_layout(L["r_reads_labels"])
    
    p_reads.hover.tooltips = [
        ("Read", "@name"), 
        ("Start", "@x0"), 
        ("End", "@x1"),
        ("Type", "@type"),
        ("Mismatches", "@mismatches"),
        ("Insertions", "@insertions")
    ]
    p_reads.hover.renderers = [L["r_reads_exons"], L["r_reads_introns"], L["r_reads_bridges"], L["r_reads_deletion_markers"]]
    p_reads.grid.grid_line_color = None
    p_reads.yaxis.visible = False
    p_reads.xaxis.formatter = NumeralTickFormatter(format="0")
    p_reads.xaxis.major_label_orientation = math.pi/4
    
    p_reads.legend.location = "top_right"
    p_reads.legend.label_text_font_size = "8pt"
    p_reads.legend.orientation = "vertical"
    
    apply_reads_visibility(p_reads, L, state)
    add_crosshair_to_plot(p_reads, L, state, plot_type="reads")
    return p_reads

def apply_reads_visibility(p_reads, L, state):
    """Syncs plot visibility settings from the UI checkboxes."""
    if "chk_xaxis_reads" in L:
        p_reads.xaxis.visible = (0 in L["chk_xaxis_reads"].active)
    if "chk_legend_reads" in L:
        p_reads.legend.visible = (0 in L["chk_legend_reads"].active)
    if "chk_show_gap_size" in L and "r_reads_bridges_labels" in L:
        L["r_reads_bridges_labels"].visible = (0 in L["chk_show_gap_size"].active)
    
    # Handle deletion markers visibility
    if "r_reads_deletion_markers" in L:
        if "chk_show_deletions" in L:
            L["r_reads_deletion_markers"].visible = (0 in L["chk_show_deletions"].active)
        else:
            # Fallback to config if checkbox not initialized yet
            show_del = state["config"]["full_read_layout"].get("show_deletion_markers", True)
            L["r_reads_deletion_markers"].visible = show_del

def apply_reads_styles(L, styles):
    """Updates glyph styles from the global configuration."""
    if "r_reads_exons" in L:
        L["r_reads_exons"].glyph.line_alpha = styles.get("exon_single", {}).get("alpha", 0.6)
    if "r_reads_bridges" in L:
        mb_s = styles.get("mate_bridge", {})
        L["r_reads_bridges"].glyph.line_color = mb_s.get("color", "red")
        L["r_reads_bridges"].glyph.line_width = mb_s.get("line_width", 1)
        L["r_reads_bridges"].glyph.line_alpha = mb_s.get("alpha", 0.4)
        mb_dash = mb_s.get("line_dash", "solid")
        if isinstance(mb_dash, list): mb_dash = tuple(mb_dash)
        L["r_reads_bridges"].glyph.line_dash = mb_dash
    if "r_reads_introns" in L:
        it_s = styles.get("intron_target", {})
        it_dash = it_s.get("line_dash", "solid")
        if isinstance(it_dash, list): it_dash = tuple(it_dash)
        L["r_reads_introns"].glyph.line_dash = it_dash
    if "r_reads_deletion_markers" in L:
        del_s = styles.get("deletion", {})
        col = del_s.get("color", "red")
        L["r_reads_deletion_markers"].glyph.fill_color = col
        L["r_reads_deletion_markers"].glyph.line_color = col
        L["r_reads_deletion_markers"].glyph.size = del_s.get("marker_size", 12)
        L["r_reads_deletion_markers"].glyph.line_alpha = del_s.get("alpha", 0.9)
        L["r_reads_deletion_markers"].glyph.line_width = 2
    if "r_reads_labels" in L:
        sl_s = styles.get("label_sample", {})
        L["r_reads_labels"].text_color = sl_s.get("color", "white")
        L["r_reads_labels"].text_font_size = sl_s.get("font_size", "10pt")
        L["r_reads_labels"].background_fill_color = sl_s.get("background_color", "black")
        L["r_reads_labels"].background_fill_alpha = sl_s.get("background_alpha", 1.0)
