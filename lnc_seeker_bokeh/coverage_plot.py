# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import math
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource, HoverTool, Range1d, WheelZoomTool, PanTool, 
    FixedTicker, NumeralTickFormatter, LabelSet, CustomJS, LogTicker
)
from bokeh.events import Reset
import lnc_seeker
from lnc_seeker_bokeh.state import log_safe, strip_id
from lnc_seeker_bokeh.plotting_base import add_crosshair_to_plot
from lnc_seeker_bokeh.data_utils import categorize_junction

def update_sample_coverage(state, res, s_data, L, name=None):
    """Updates only the coverage data source for a plot, with automatic resampling based on view range."""
    import numpy as np
    import lnc_seeker
    from lnc_seeker_bokeh.state import strip_id

    # res = (p, ds_cov)
    p, ds_cov = res[0], res[1]

    positions = s_data.get("positions", [])
    if len(positions) == 0:
        ds_cov.data = dict(x=[], y_bg=[], y_fg=[])
        return

    depths = np.array(s_data.get("depths", []))
    depths_hq = np.array(s_data.get("depths_hq", []))

    show_full_cov = (0 in L["chk_full_cov"].active) if "chk_full_cov" in L else False
    show_bg = (0 in L["chk_show_bg"].active) if "chk_show_bg" in L else False
    normalize = (0 in L["chk_normalize"].active) if "chk_normalize" in L else False

    factor = 1.0
    if normalize and name:
        raw_name = strip_id(name)
        gene = L["sel_gene"].value if "sel_gene" in L else None
        if gene and gene in state.get("bam_hierarchy", {}):
            cohorts = state["bam_hierarchy"][gene].get("cohorts", {})
            if raw_name in cohorts:
                factor = float(cohorts[raw_name].get("num_samples", 1.0))
                if factor <= 0: factor = 1.0

    if not show_full_cov and not show_bg:
        ds_cov.data = dict(x=[], y_bg=[], y_fg=[])
        return

    # 1. Determine Target Resolution (Pixel Width of the diagram)
    fixed_width = (0 in L["chk_fixed_width"].active) if L.get("chk_fixed_width") else False
    if fixed_width and L.get("spn_plot_width"):
        pixel_width = int(L["spn_plot_width"].value)
    else:
        # p.inner_width is standard for Bokeh figures; fallback if not yet rendered
        pixel_width = getattr(p, 'inner_width', 0)
        if not pixel_width or pixel_width <= 0:
            pixel_width = 1200

    # 2. Get the current viewport range
    try:
        cur_start, cur_end = p.x_range.start, p.x_range.end
        if cur_start is None or cur_end is None or cur_start >= cur_end:
             cur_start, cur_end = positions[0], positions[-1]
    except Exception:
        cur_start, cur_end = positions[0], positions[-1]

    # 3. Filter data to the visible domain + small margin
    idx_start = np.searchsorted(positions, cur_start)
    idx_end = np.searchsorted(positions, cur_end)

    margin = 50
    idx_start = max(0, idx_start - margin)
    idx_end = min(len(positions), idx_end + margin)

    vis_pos = positions[idx_start:idx_end]
    vis_bg = depths[idx_start:idx_end]
    vis_fg = depths_hq[idx_start:idx_end]

    y_bg_full = np.maximum(vis_bg / factor, 0.1)
    y_fg_full = np.maximum(vis_fg / factor, 0.1)

    # 4. Downsample to match pixel resolution (1 point per pixel)
    if len(vis_pos) <= pixel_width:
        ds_cov.data = dict(x=vis_pos, y_bg=y_bg_full, y_fg=y_fg_full)
    else:
        try:
            all_pts = s_data.get("junction_points", [])
            cliff_positions = [int(pt['position']) for pt in all_pts if cur_start <= pt['position'] <= cur_end]

            x_list, ybg_list, yfg_list = lnc_seeker.downsample_coverage_py(
                list(vis_pos), y_bg_full.tolist(), y_fg_full.tolist(),
                int(pixel_width), cliff_positions
            )
            ds_cov.data = dict(x=x_list, y_bg=ybg_list, y_fg=yfg_list)
        except Exception as e:
            log_safe(state, f"Auto-resampling failed for {res[0].title.text}: {e}")
            ds_cov.data = dict(x=vis_pos, y_bg=y_bg_full, y_fg=y_fg_full)

def update_sample_data(state, name, res, s_data, L, min_reads, active_types, show_marked, marked_sets):
    """Updates the data sources for a coverage plot."""
    import numpy as np
    from lnc_seeker_bokeh.constants import RedGrayBlue11
    from lnc_seeker_bokeh.plotting_base import get_j_color
    from lnc_seeker_bokeh.state import strip_id

    p, ds_cov, ds_js, ds_jp, r_bg, r_fg, r_arcs, r_cliffs, ds_cliffs, r_cliff_steps, ds_cliff_steps, ds_cliff_area, r_cliff_area, r_ambiguity, ds_ambiguity, handler = res
    
    # Normalization factor
    normalize = (0 in L["chk_normalize"].active) if "chk_normalize" in L else False
    normalize_js = (0 in L["chk_normalize_junctions"].active) if "chk_normalize_junctions" in L else True
    
    factor = 1.0
    if normalize:
        raw_name = strip_id(name)
        gene = L["sel_gene"].value if "sel_gene" in L else None
        if gene and gene in state.get("bam_hierarchy", {}):
            cohorts = state["bam_hierarchy"][gene].get("cohorts", {})
            if raw_name in cohorts:
                factor = float(cohorts[raw_name].get("num_samples", 1.0))
                if factor <= 0: factor = 1.0

    j_factor = factor if normalize_js else 1.0
    
    positions = s_data.get("positions", [])
    depths = np.array(s_data.get("depths", []))
    depths_hq = np.array(s_data.get("depths_hq", []))
    depths_amb = np.array(s_data.get("depths_ambiguity", depths_hq))
    
    show_full_cov = (0 in L["chk_full_cov"].active) if "chk_full_cov" in L else False
    show_bg = (0 in L["chk_show_bg"].active) if "chk_show_bg" in L else False

    if len(positions) > 0:
        try:
            update_sample_coverage(state, res, s_data, L, name=name)
        except Exception as e: log_safe(state, f"Error preparing coverage data for {name}: {e}")

        max_depth = float(np.max(depths)) / factor
        amb_x0, amb_x1, amb_y0, amb_y1 = [], [], [], []
        amb_cfg = state["config"]["coverage_and_junctions_profile"].get("high_ambiguity_highlighting", {}).get("ambiguity_highlight", {})
        thresh = amb_cfg.get("threshold", 3.0)
        is_ambig = (depths > depths_amb * thresh) & (depths > 2) 
        if np.any(is_ambig):
            diff = np.diff(is_ambig.astype(int), prepend=0, append=0)
            starts = np.where(diff == 1)[0]
            ends = np.where(diff == -1)[0] - 1
            for s, e in zip(starts, ends):
                e_idx = min(e + 1, len(positions) - 1)
                amb_x0.append(positions[s]); amb_x1.append(positions[e_idx] if e_idx > s else positions[s] + 1)
                amb_y0.append(0.01); amb_y1.append(1000000)
        ds_ambiguity.data = dict(x0=amb_x0, x1=amb_x1, y0=amb_y0, y1=amb_y1)

    filter_clean = (0 in L["chk_filter_flanks"].active) if "chk_filter_flanks" in L else False
    spans = s_data.get("junction_spans", [])
    x0, y0, x1, y1, cx, cy, mid, h_peak, label_w, label_h, marker_line, colors, labels, line_ws, refs, fsizes, fstyles, anchored_start_list, anchored_end_list = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    for row_sj in spans:
        try:
            reads = row_sj.get('reads_clean', row_sj['reads']) if filter_clean else row_sj['reads']
            if reads < min_reads: continue

            start, end = row_sj['start'], row_sj['end']
            j = (start, end)
            j_type = categorize_junction(state, start, end)
            is_marked = j in marked_sets.get(j_type, set())
            if is_marked:
                if not show_marked: continue
                color_key = "marked"
            else:
                if j_type not in active_types: continue
                color_key = j_type
            x0.append(start); x1.append(end); y0.append(0.1); y1.append(0.1)
            anchored_start_list.append(row_sj.get('anchored_start', 0)); anchored_end_list.append(row_sj.get('anchored_end', 0))
            mid_val = (start+end)/2; mid.append(mid_val); cx.append(mid_val)
            cy_val = (float(reads) / j_factor) * 1.5 + 2.0; cy.append(cy_val)
            
            h_peak.append(math.sqrt(cy_val * 0.1))
            labels.append(f"{reads/j_factor:.1f}" if j_factor > 1.0 else str(int(reads)))
            
            label_w.append(20 * (2 if reads/j_factor > 99 else 1)); label_h.append(14); colors.append(get_j_color(state, color_key, 0)); marker_line.append(get_j_color(state, color_key, 1))
            line_ws.append(max(1, min(8, int(math.log10(reads/j_factor + 1) * 2) + 1))); refs.append(row_sj.get('reference', '')); fsizes.append("8pt"); fstyles.append("normal")
        except Exception as e:
            log_safe(state, f"Error processing junction span: {e}")
            continue
    ds_js.data = dict(x0=x0, y0=y0, x1=x1, y1=y1, cx=cx, cy=cy, mid=mid, h_peak=h_peak, label_w=label_w, label_h=label_h, marker_line=marker_line, color=colors, line_w=line_ws, label=labels, reference=refs, fsize=fsizes, fstyle=fstyles, anchored_start=anchored_start_list, anchored_end=anchored_end_list)

    pts = s_data.get("junction_points", [])
    px, py, p_colors, p_alphas, p_labels = [], [], [], [], []
    cliff_x, cliff_y0, cliff_y1, cliff_colors, cliff_labels = [], [], [], [], []
    for row_p in pts:
        change = row_p.get('change_pct')
        if not math.isfinite(change): continue
        pos = row_p['position']
        avg_before = max(row_p['avg_before'] / factor, 0.1)
        avg_after = max(row_p['avg_after'] / factor, 0.1)
        px.append(pos)
        py.append(avg_before)
        idx = int((change + 100) / 200 * (len(RedGrayBlue11)-1))
        idx = max(0, min(len(RedGrayBlue11)-1, idx))
        p_colors.append(RedGrayBlue11[idx])
        p_alphas.append(max(0.5, min(1.0, abs(change)/100.0)))
        p_labels.append(change)
        cliff_x.append(pos)
        cliff_y0.append(avg_before)
        cliff_y1.append(avg_after)
        cliff_colors.append(RedGrayBlue11[idx])
        cliff_labels.append(change)
    ds_jp.data = dict(x=px, y=py, color=p_colors, alpha=p_alphas, label=p_labels)
    ds_cliffs.data = dict(x=cliff_x, y0=cliff_y0, y1=cliff_y1, color=cliff_colors, label=cliff_labels)
    
    h_x0, h_x1, h_y0, h_y1 = [], [], [], []
    sorted_sig = sorted(pts, key=lambda x: x['position'])
    for i in range(len(sorted_sig)-1):
        p1, p2 = sorted_sig[i], sorted_sig[i+1]
        y_p1_after = max(p1['avg_after'] / factor, 0.1)
        y_p2_before = max(p2['avg_before'] / factor, 0.1)
        mid_x = (p1['position'] + p2['position']) / 2
        h_x0.append(p1['position']); h_x1.append(mid_x); h_y0.append(y_p1_after); h_y1.append(y_p1_after)
        h_x0.append(mid_x); h_x1.append(mid_x); h_y0.append(y_p1_after); h_y1.append(y_p2_before)
        h_x0.append(mid_x); h_x1.append(p2['position']); h_y0.append(y_p2_before); h_y1.append(y_p2_before)
    ds_cliff_steps.data = dict(x0=h_x0, x1=h_x1, y0=h_y0, y1=h_y1)
    
    area_x, area_y = [], []
    if sorted_sig:
        for i in range(len(sorted_sig)):
            pt = sorted_sig[i]; y_before = max(pt['avg_before'] / factor, 0.1); y_after = max(pt['avg_after'] / factor, 0.1); pos = pt['position']
            area_x.extend([pos, pos]); area_y.extend([y_before, y_after])
            if i < len(sorted_sig) - 1:
                p_next_sig = sorted_sig[i+1]; mid_x_sig = (pos + p_next_sig['position']) / 2; y_next_before_sig = max(p_next_sig['avg_before'] / factor, 0.1)
                area_x.append(mid_x_sig); area_y.append(y_after); area_x.append(mid_x_sig); area_y.append(y_next_before_sig); area_x.append(p_next_sig['position']); area_y.append(y_next_before_sig)
    ds_cliff_area.data = dict(x=area_x, y=area_y)

def create_sample_plot(state, sample_name, L, on_reset_zoom_click, sync_all_junction_selections, x_range=None):

    """Creates a coverage/junction plot for a sample."""
    target_h = L["spn_height"].value if "spn_height" in L else state["config"]["coverage_and_junctions_profile"].get("plot_height", 450)
    
    # Extract metadata for the title
    tissue, status, num_samples = L["selection_manager"].get_sample_metadata(sample_name)

    p = figure(
        height=target_h, 
        title=f"Coverage & Junctions Profile - {sample_name} ({tissue}, {status}, samples={num_samples})",
        x_axis_label="Genomic Position", y_axis_label="Depth (Log Scale)",
        y_axis_type="log",
        sizing_mode="stretch_width",
        tools="reset,save,tap",
        min_border_left=80,
        x_range=x_range if x_range is not None else Range1d(start=0, end=1),
        output_backend=state["config"]["general"].get("output_backend", "canvas")
    )
    # Lock tools to x-axis only
    w_zoom = WheelZoomTool(dimensions="width")
    p.add_tools(w_zoom, PanTool(dimensions="width"))
    p.toolbar.active_scroll = w_zoom
    p.on_event(Reset, on_reset_zoom_click)
    p.js_on_event(Reset, CustomJS(args=dict(xr=p.x_range, yr=p.y_range, dr=L["ds_core_range"]), code="""
        if (dr.data['start'] && dr.data['start'].length > 0) {
            xr.start = dr.data['start'][0];
            xr.end = dr.data['end'][0];
        }
        if (dr.data['y_max'] && dr.data['y_max'].length > 0) {
            yr.start = 0.1;
            yr.end = dr.data['y_max'][0];
        }
    """))

    # ... (Wait, I'll just find the end of the create_sample_plot function)
    
    # Sources
    ds_cov = ColumnDataSource(data=dict(x=[], y_bg=[], y_fg=[]))
    ds_cliffs = ColumnDataSource(data=dict(x=[], y0=[], y1=[], color=[], label=[]))
    ds_cliff_steps = ColumnDataSource(data=dict(x0=[], x1=[], y0=[], y1=[]))
    ds_cliff_area = ColumnDataSource(data=dict(x=[], y=[]))
    ds_ambiguity = ColumnDataSource(data=dict(x0=[], x1=[], y0=[], y1=[]))
    ds_j_spans = ColumnDataSource(data=dict(
        x0=[], y0=[], x1=[], y1=[], cx=[], cy=[], mid=[], h_peak=[], 
        label_w=[], label_h=[], marker_line=[], color=[], line_w=[], 
        label=[], fsize=[], fstyle=[]
    ))
    ds_j_pts = ColumnDataSource(data=dict(x=[], y=[], color=[], alpha=[], label=[]))
    
    # Coverage Tracks
    cov_styles = state["config"]["coverage_and_junctions_profile"]["coverage_components"].get("coverage_styles", {})
    bg_style = cov_styles.get("background", {"color": "#444444", "line_width": 1, "alpha": 0.5, "label": "All Reads (BG)"})
    fg_style = cov_styles.get("foreground", {"color": "navy", "line_width": 1.5, "alpha": 0.95, "label": "HQ Reads (FG)"})

    r_bg = p.step(x='x', y='y_bg', source=ds_cov, 
                  color=bg_style.get("color", "#444444"), 
                  line_width=bg_style.get("line_width", 1), 
                  alpha=bg_style.get("alpha", 0.5), 
                  mode="center", 
                  legend_label=bg_style.get("label", "All Reads (BG)"))
    
    r_fg = p.step(x='x', y='y_fg', source=ds_cov, 
                  color=fg_style.get("color", "navy"), 
                  line_width=fg_style.get("line_width", 1.5), 
                  alpha=fg_style.get("alpha", 0.95), 
                  mode="center", 
                  legend_label=fg_style.get("label", "HQ Reads (FG)"))
    
    # Ambiguity Highlight
    amb_cfg = state["config"]["coverage_and_junctions_profile"]["high_ambiguity_highlighting"].get("ambiguity_highlight", {})
    r_ambiguity = p.quad(top=1000000, bottom=0.01, left='x0', right='x1', source=ds_ambiguity,
                         fill_color=amb_cfg.get("color", "yellow"),
                         fill_alpha=amb_cfg.get("alpha", 0.1),
                         line_color=None, visible=False, legend_label="High Ambiguity")
    
    # Cliff Track (Simplified Coverage)
    opt_cfg = state["config"]["coverage_and_junctions_profile"]["coverage_components"].get("optimized_coverage", {})
    r_cliffs = p.segment(x0='x', y0='y0', x1='x', y1='y1', source=ds_cliffs, 
                         color='color', 
                         line_width=opt_cfg.get("cliff_width", 3), 
                         alpha=opt_cfg.get("cliff_alpha", 0.8), 
                         legend_label="Coverage Cliffs", visible=False)
    
    r_cliff_area = p.varea(x='x', y1='y', y2=0.1, source=ds_cliff_area,
                           fill_color=opt_cfg.get("line_color", "gray"),
                           fill_alpha=opt_cfg.get("area_alpha", 0.1),
                           legend_label="Cliff Trend", visible=False)

    r_cliff_steps = p.segment(x0='x0', y0='y0', x1='x1', y1='y1', source=ds_cliff_steps, 
                              color=opt_cfg.get("line_color", "gray"), 
                              line_width=opt_cfg.get("line_width", 1), 
                              alpha=opt_cfg.get("trend_alpha", 0.4), 
                              legend_label="Cliff Trend", visible=False)
    
    # Junction Arcs
    r_arcs = p.quadratic(x0='x0', y0='y0', x1='x1', y1='y1', cx='cx', cy='cy', 
                source=ds_j_spans, line_width='line_w', color='color', alpha=0.7, legend_label="Splice Junctions")
    r_arcs.selection_glyph = None
    r_arcs.nonselection_glyph = r_arcs.glyph.clone(line_alpha=0.6)
    
    # Peak Marker (Square)
    marker_renderer = p.scatter(x='mid', y='h_peak', size=12, fill_color='color', line_color='marker_line', alpha=0.95, marker="square", source=ds_j_spans, legend_label="Junction Read")
    marker_renderer.selection_glyph = None
    marker_renderer.nonselection_glyph = marker_renderer.glyph.clone(fill_alpha=0.9, line_alpha=0.9)

    def on_junction_select(attr, old, new):
        if L["selection_updating"]:
            return
            
        if old != new:
            L["is_squeezing"] = False
            
        # Update global pointers so sliders use the correct plot context
        L["ds_j_spans"] = ds_j_spans
        L["on_junction_select"] = on_junction_select

        if not new:
            # Reset all to default if selection cleared
            L["last_selected_junctions"] = []
            L["last_junction"] = None
            L["selection_updating"] = True
            try:
                for s_name, plot_objs in L["sample_plots"].items():
                    ds_js = plot_objs[2]
                    if ds_js.selected.indices:
                        ds_js.selected.indices = []
                    
                    data = ds_js.data
                    n = len(data['label'])
                    data['label_w'] = [20 * (2 if float(r) > 99 else 1) for r in data['label']]
                    data['label_h'] = [14] * n
                    data['fsize'] = ["8pt"] * n
                    data['fstyle'] = ["normal"] * n
                    data['line_w'] = [max(1, min(8, int(math.log10(float(r) + 1) * 2) + 1)) for r in data['label']]
                    ds_js.trigger('data', ds_js.data, ds_js.data)
                
                # Clear Transcript Creator elements too
                L["mul_curr_j"].value = []
                if L.get("ds_creator_annotations"):
                    L["ds_creator_annotations"].selected.indices = []
            finally:
                L["selection_updating"] = False
            return

        data = dict(ds_j_spans.data)
        selected_junctions = []
        for idx in new:
            selected_junctions.append((data['reference'][idx], data['x0'][idx], data['x1'][idx]))
        
        if not selected_junctions:
            return

        # Store selected junctions for restoration/refresh
        L["last_selected_junctions"] = list(selected_junctions)
        ref_last, start_last, end_last = selected_junctions[-1]
        L["last_junction"] = (ref_last, start_last, end_last)

        # Master synchronization across all components
        sync_all_junction_selections(source="sample_plot", new=selected_junctions)
        
        target_samples_bams = []
        for s_name in L["sel_samples"].value:
            bp = L["selection_manager"].get_bam_path(s_name)
            if bp:
                target_samples_bams.append((s_name, bp))
        L["last_target_samples_bams"] = target_samples_bams

        # Use unified collective fetch from ReadsManager
        if "reads_manager" in L:
            L["reads_manager"].update_reads_ui(selected_junctions, target_samples_bams)

    ds_j_spans.selected.on_change('indices', on_junction_select)
    L["ds_j_spans"] = ds_j_spans
    L["on_junction_select"] = on_junction_select

    # Read Count Labels
    r_labels = p.rect(x='mid', y='h_peak', width='label_w', height='label_h', width_units='screen', height_units='screen', source=ds_j_spans, fill_color='marker_line', fill_alpha=0.95, line_color=None)
    r_labels.selection_glyph = None
    r_labels.nonselection_glyph = r_labels.glyph.clone(fill_alpha=0.9)
    
    p.add_layout(LabelSet(x='mid', y='h_peak', text='label', source=ds_j_spans, text_font_size="fsize", text_font_style="fstyle", text_color="white", text_align="center", text_baseline="middle"))
    
    # Junction Points (Change)
    change_renderer = p.scatter(x='x', y='y', source=ds_j_pts, size=10, fill_color='color', line_color="black", line_width=1, alpha='alpha', legend_label="Junction Change")
    p.add_tools(HoverTool(renderers=[change_renderer], tooltips=[('Pos', '@x'), ('Change', '@label{0.0}%')], mode='mouse'))

    hover_arc = HoverTool(renderers=[marker_renderer, r_labels], tooltips=[
        ("Reference", "@x0 to @x1"),
        ("Spanning Reads", "@label"),
        ("Anchored (5')", "@anchored_start"),
        ("Anchored (3')", "@anchored_end")
    ], mode='mouse')
    p.add_tools(hover_arc)
    
    p.legend.location = "top_right"
    p.legend.orientation = "vertical"
    p.legend.click_policy = "hide"
    p.legend.label_text_font_size = "8pt"
    
    hover_cov = HoverTool(renderers=[r_fg], tooltips=[("Pos", "@x"), ("Depth", "@y_fg")])
    p.add_tools(hover_cov)

    # Add Hover for Cliffs
    hover_cliffs = HoverTool(renderers=[r_cliffs], tooltips=[("Pos", "@x"), ("Change", "@label{0.0}%")])
    p.add_tools(hover_cliffs)

    p.xaxis.formatter = NumeralTickFormatter(format="0")
    p.yaxis.formatter = NumeralTickFormatter(format="0")
    p.xaxis.major_label_orientation = math.pi/4
    
    # Gene-region extension boundaries if present (dictionary-based selection)
    if "ds_extension" in L:
        print(f"DEBUG: Adding extension boundaries to sample plot for {sample_name}")
        r_ext = p.segment(x0='x0', y0='y0', x1='x1', y1='y1', source=L["ds_extension"],
                          line_color='color', line_width='width', line_dash='dash', alpha='alpha',
                          name="ext_segment")
        p.add_tools(HoverTool(renderers=[r_ext], tooltips=[
            ("Boundary", "@desc"),
            ("Extraction Region", "@region"),
            ("Position", "@x0{0,0}")
        ]))

    add_crosshair_to_plot(p, L, state, plot_type="coverage")
    return (p, ds_cov, ds_j_spans, ds_j_pts, r_bg, r_fg, r_arcs, r_cliffs, ds_cliffs, r_cliff_steps, ds_cliff_steps, ds_cliff_area, r_cliff_area, r_ambiguity, ds_ambiguity, on_junction_select)

def get_or_update_sample_plot(state, name, L, on_reset, on_selection, 
                              cov_target_h, backend, fixed_width, target_width, 
                              locked_y_top, show_full_cov, show_bg, show_cliffs):
    is_new_plot = False
    if name not in L["sample_plots"]:
        res = create_sample_plot(state, name, L, on_reset, on_selection, x_range=L["shared_x_range"])
        p, ds_cov, ds_js, ds_jp, r_bg, r_fg, r_arcs, r_cliffs, ds_cliffs, r_cliff_steps, ds_cliff_steps, ds_cliff_area, r_cliff_area, r_ambiguity, ds_ambiguity, handler = res
        is_new_plot = True
        p.y_range = Range1d(start=0.1, end=locked_y_top)
        p.height = cov_target_h
        p.output_backend = backend
        if fixed_width: p.sizing_mode = "fixed"; p.width = target_width
        else: p.sizing_mode = "stretch_width"
        L["sample_plots"][name] = res
        
        p.on_event(Reset, on_reset)
        p.js_on_event(Reset, CustomJS(args=dict(xr=p.x_range, yr=p.y_range, dr=L["ds_core_range"]), code="""
            if (dr.data['start'] && dr.data['start'].length > 0) { xr.start = dr.data['start'][0]; xr.end = dr.data['end'][0]; }
            if (dr.data['y_max'] && dr.data['y_max'].length > 0) { yr.start = 0.1; yr.end = dr.data['y_max'][0]; }
        """))
    else:
        res = L["sample_plots"][name]
        p, ds_cov, ds_js, ds_jp, r_bg, r_fg, r_arcs, r_cliffs, ds_cliffs, r_cliff_steps, ds_cliff_steps, ds_cliff_area, r_cliff_area, r_ambiguity, ds_ambiguity, handler = res
        
        # Ensure extension boundaries are added to existing plots (migration for live reload)
        has_ext = any(r.name == 'ext_segment' for r in p.renderers)
        if not has_ext and "ds_extension" in L:
            r_ext = p.segment(x0='x0', y0='y0', x1='x1', y1='y1', source=L["ds_extension"],
                              line_color='color', line_width='width', line_dash='dash', alpha='alpha',
                              name="ext_segment")
            p.add_tools(HoverTool(renderers=[r_ext], tooltips=[
                ("Boundary", "@desc"), ("Extraction Region", "@region"), ("Position", "@x0{0,0}")
            ]))

        if p.x_range != L["shared_x_range"]: p.x_range = L["shared_x_range"]
        if abs(p.y_range.end - locked_y_top) > 0.1:
            p.y_range.start = 0.1; p.y_range.end = locked_y_top
        p.height = cov_target_h
        p.output_backend = backend
        if fixed_width: p.sizing_mode = "fixed"; p.width = target_width
        else: p.sizing_mode = "stretch_width"
    
    r_fg.visible = show_full_cov
    r_bg.visible = show_bg
    r_cliffs.visible = show_cliffs
    r_cliff_steps.visible = show_cliffs
    r_cliff_area.visible = show_cliffs
    r_ambiguity.visible = (0 in L["chk_ambiguity"].active) if "chk_ambiguity" in L else False
    
    if "chk_xaxis_main" in L: p.xaxis.visible = (0 in L["chk_xaxis_main"].active)
    if "chk_legend_main" in L: p.legend.visible = (0 in L["chk_legend_main"].active)

    # Update Y-Axis Styles based on normalization
    normalize = (0 in L["chk_normalize"].active) if "chk_normalize" in L else False
    y_styles = state.get("config", {}).get("coverage_and_junctions_profile", {}).get("y_axis_styles", {})
    style_key = "normalized" if normalize else "absolute"
    style = y_styles.get(style_key, {})
    
    if style:
        sub_text = style.get("sub_text", style_key)
        p.yaxis.axis_label = f"Depth {sub_text}"
        p.yaxis.axis_label_text_align = "center"
        
        if "label_text_color" in style: p.yaxis.axis_label_text_color = style["label_text_color"]
        if "major_label_text_color" in style: p.yaxis.major_label_text_color = style["major_label_text_color"]
        if "major_tick_line_color" in style: p.yaxis.major_tick_line_color = style["major_tick_line_color"]
        if "minor_tick_line_color" in style: p.yaxis.minor_tick_line_color = style["minor_tick_line_color"]
    
    # Update tickers based on normalization
    if normalize:
        p.yaxis.ticker = LogTicker()
    else:
        # 1 Mio = 1,000,000
        p.yaxis.ticker = FixedTicker(ticks=[1, 20, 50, 100, 1000, 10000, 100000, 1000000])
        
    return res, is_new_plot
