# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
from bokeh.models import Span, CustomJS

def get_j_color(state, j_type, variant=0):
    jc = state.get("config", {}).get("junctions_splicing", {}).get("junction_colors", {})
    # Try exact key, then normalized variants (spaces <-> underscores), then lowercase
    candidates = [
        j_type,
        j_type.replace(" ", "_"),
        j_type.replace("_", " "),
        j_type.lower(),
        j_type.lower().replace(" ", "_")
    ]
    colors = None
    for k in candidates:
        if k in jc:
            colors = jc[k]
            break
    if not colors:
        colors = ["gray", "black"]
        
    return colors[variant] if variant < len(colors) else colors[0]

def get_transcript_color(state, transcript_id, feature="exon"):
    tc = state.get("config", {}).get("genome_annotations", {}).get("transcript_colors", {})
    if feature == "intron":
        return tc.get("intron", "black")
    
    # Check for prefix matches in the config (e.g., NM_, NR_, XM_, XR_)
    if transcript_id and isinstance(transcript_id, str):
        # Sort prefixes by length descending to ensure the most specific match first
        prefixes = sorted([k for k in tc.keys() if k not in ("default", "intron")], key=len, reverse=True)
        for prefix in prefixes:
            if transcript_id.startswith(prefix):
                return tc[prefix]
                
    return tc.get("default", "steelblue")

def get_status_style(state, status):
    """Returns a style dictionary for the status badge based on config."""
    sc = state.get("config", {}).get("general", {}).get("status_colors", {})
    
    # Normalize status for lookup
    s_key = status.lower() if status else "default"
    
    # Try exact match, then default
    style = sc.get(s_key)
    if not style:
        style = sc.get("default", { "bg": "#e1f5fe", "fg": "#01579b" })
        
    return style

def add_crosshair_to_plot(p, L, state, plot_type="unknown"):
    """Adds a crosshair (vertical and horizontal lines) that is synchronized across all plots."""
    # Safety check: avoid adding multiple crosshair spans to the same plot
    for r in p.renderers:
        if isinstance(r, Span) and getattr(r, "name", "") in ["crosshair_v", "crosshair_h", "crosshair_span"]:
            return

    c_cfg = state.get("config", {}).get("general", {}).get("crosshair", {})
    c_color = c_cfg.get("color", "red")
    c_color_h = c_cfg.get("horizontal_color", c_color)
    c_width = c_cfg.get("width", 1)
    c_alpha = c_cfg.get("alpha", 0.8)
    
    # Vertical line
    span_v = Span(dimension="height", line_dash="solid", 
                  line_width=c_width, line_color=c_color, line_alpha=c_alpha,
                  location=0, visible=False, level="overlay", name="crosshair_v")
    
    # Horizontal line
    span_h = Span(dimension="width", line_dash="solid", 
                  line_width=c_width, line_color=c_color_h, line_alpha=c_alpha,
                  location=0, visible=False, level="overlay", name="crosshair_h")
    
    p.add_layout(span_v)
    p.add_layout(span_h)
    
    L["ds_crosshair"].js_on_change("data", CustomJS(args=dict(p=p, ds=L["ds_crosshair"], span_v=span_v, span_h=span_h, chk=L["chk_crosshair"], chk_ext=L["chk_crosshair_ext"], current_type=plot_type), code="""
        if (ds.data["x"] && ds.data["x"].length > 0 && ds.data["y"] && ds.data["y"].length > 0) {
            const active = chk.active.includes(0);
            const extended = chk_ext.active.includes(0);
            const source_type = ds.data["type"][0];
            const source_id = ds.data["id"] ? ds.data["id"][0] : null;
            
            span_v.location = ds.data["x"][0];
            span_h.location = ds.data["y"][0];
            
            const visible = active && ds.data["visible"][0];
            span_v.visible = visible;
            
            // Standard mode: horizontal line only in source plot
            // Extended mode: horizontal line in source plot OR all coverage plots if source is coverage
            if (extended && source_type === "coverage" && current_type === "coverage") {
                span_h.visible = visible;
            } else {
                // Only show horizontal line if this is the EXACT plot instance the mouse is in
                span_h.visible = visible && (source_id === p.id);
            }
        }
    """))

    p.js_on_event("mousemove", CustomJS(args=dict(p=p, ds=L["ds_crosshair"], chk=L["chk_crosshair"], plot_type=plot_type), code="""
        if (cb_obj.x !== undefined && cb_obj.x !== null && cb_obj.y !== undefined && cb_obj.y !== null) {
            const active = chk.active.includes(0);
            if (active) {
                ds.data = { x: [cb_obj.x], y: [cb_obj.y], type: [plot_type], id: [p.id], visible: [true] };
            }
        }
    """))
