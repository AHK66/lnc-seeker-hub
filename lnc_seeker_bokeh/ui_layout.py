# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
from bokeh.models import (
    ColumnDataSource, Range1d, Div, Select, MultiSelect, Spinner, 
    CheckboxGroup, Slider, TextInput, NumericInput, Button,
    Tabs, TabPanel, FileInput
)
from bokeh.layouts import column, row
from lnc_seeker_bokeh.state import log_safe
import importlib.resources as resources
from pathlib import Path

def initialize_local_state(app):
    state = app.state
    config = state.get("config")
    if config is None:
        # Fallback to prevent crash if config.json failed to load
        config = {"data_selection": {}, "genome_annotations": {}, "coverage_and_junctions_profile": {}, "full_read_layout": {}, "transcript_creator": {}, "general": {}}
    
    genes = sorted(list(state.get("bam_hierarchy", {}).keys()))
    
    # Get GTF options
    gtf_paths = config.get("data_selection", {}).get("gtf_paths", [])
    if "selected_gtfs" not in config.get("data_selection", {}):
        if "data_selection" not in config: config["data_selection"] = {}
        config["data_selection"]["selected_gtfs"] = [gtf_paths[0]] if gtf_paths else []
    
    selected_gtfs = config["data_selection"]["selected_gtfs"]
    
    # Ensure the first GTF is always selected
    if gtf_paths:
        if gtf_paths[0] not in selected_gtfs:
            selected_gtfs.insert(0, gtf_paths[0])
            config["data_selection"]["selected_gtfs"] = selected_gtfs

    active_gtf_indices = [i for i, p in enumerate(gtf_paths) if p in selected_gtfs]
    if not active_gtf_indices and gtf_paths:
        active_gtf_indices = [0]

    # Load external manual HTML (prefer importlib.resources; fallback to local file during development)
    try:
        try:
            _manual_html = resources.files("lnc_seeker_bokeh").joinpath("static", "manual.html").read_text(encoding="utf-8")
        except Exception:
            _manual_html = resources.read_text("lnc_seeker_bokeh", "static/manual.html")
    except Exception:
        _manual_html = (Path(__file__).resolve().parent / "static" / "manual.html").read_text(encoding="utf-8")

    L = {
        "ds_annotations": ColumnDataSource(data=dict(left=[], right=[], mid=[], width=[], y=[], label=[], transcript=[], exon_num=[], cds_num=[], feature=[], strand=[], color=[])),
        "ds_cds": ColumnDataSource(data=dict(left=[], right=[], mid=[], width=[], y=[], label=[], transcript=[], exon_num=[], cds_num=[], feature=[], strand=[], color=[])),
        "ds_transcripts": ColumnDataSource(data=dict(start=[], end=[], y=[], strand=[], transcript=[])),
        "ds_intron_markers": ColumnDataSource(data=dict(x=[], y=[], angle=[], transcript=[], strand=[])),
        "ds_unsupported_introns": ColumnDataSource(data=dict(x=[], y=[])),
        "ds_introns": ColumnDataSource(data=dict(left=[], right=[], mid=[], width=[], y=[], transcript=[], gene_id=[], strand=[], length=[], reference=[])),
        "ds_gene_labels": ColumnDataSource(data=dict(x=[], y=[], text=[])),
        "ds_extension": ColumnDataSource(data=dict(x=[], y0=[], y1=[], desc=[])),
        "ds_reads_exons": ColumnDataSource(data=dict(x0=[], x1=[], y=[], color=[], thickness=[], name=[], mismatches=[], insertions=[], type=[], label=[], mid_x=[], label_color=[])),
        "ds_reads_introns": ColumnDataSource(data=dict(x0=[], x1=[], y=[], label=[], mid_x=[], color=[], line_w=[], name=[])),
        "ds_reads_bridges": ColumnDataSource(data=dict(x0=[], x1=[], y=[], label=[], mid_x=[], color=[], name=[])),
        "ds_reads_full": ColumnDataSource(data=dict(x0=[], x1=[], y=[], name=[], color=[])),
        "ds_reads_labels": ColumnDataSource(data=dict(x=[], y=[], text=[], color=[])),
        "ds_crosshair": ColumnDataSource(data=dict(x=[], y=[], type=[], id=[], visible=[])),
        "ds_creator_annotations": ColumnDataSource(data=dict(left=[], right=[], mid=[], width=[], h=[], y=[], label=[], transcript=[], feature_num=[], feature=[], strand=[], color=[])),
        "ds_creator_transcripts": ColumnDataSource(data=dict(start=[], end=[], y=[], strand=[], transcript=[])),
        "ds_creator_markers": ColumnDataSource(data=dict(x=[], y=[], angle=[], transcript=[], strand=[])),
        "sample_plots": {}, # sample_name -> (p, ds_cov, ds_j_spans, ds_j_pts, ...)
        "p_ann": None,
        "p_creator": None,
        "p_reads": None,
        "shared_x_range": Range1d(start=0, end=1),
        "ds_core_range": ColumnDataSource(data=dict(start=[0], end=[1], y_max=[10])),
        "data_rendered": False,
        "update_pending": False,
        "reset_x_range": True,
        "reset_y_range_ann": True,
        "tabs": None,
        "div_title": Div(text="""
            <h2 style='margin:0'>lncRNA Seeker Hub</h2>
            <script>
                if (window.location.search.includes('bokeh-session-id')) {
                    const url = new URL(window.location);
                    url.searchParams.delete('bokeh-session-id');
                    window.history.replaceState({}, '', url);
                    console.log("[SESSION] Stripped session ID from URL to avoid reload loops");
                }
            </script>
        """, styles={'padding-bottom': '10px'}),
        "div_progress": Div(text="", sizing_mode="stretch_width", height=75),
        "div_manual": Div(text=_manual_html, sizing_mode="stretch_width"),

        "div_config_warnings": Div(text="", sizing_mode="stretch_width", visible=False),
        "div_fatal_error": Div(text="", sizing_mode="stretch_width", visible=False),
        "div_gene_info": Div(text="", sizing_mode="stretch_width", visible=False),

        "sel_gene": Select(title="Gene Selection", options=[""] + genes, value=""),
        "btn_show_cohort_selection": Button(label="Show Cohort Selection", button_type="default", disabled=True, sizing_mode="stretch_width"),
        "sel_samples": MultiSelect(title="Sample Selection (Hold Ctrl/Cmd)", options=[], value=[], size=8),
        "div_gtf_header": Div(text="<b>GTF-Selection:</b>", styles={"margin-top": "10px", "margin-bottom": "0px"}),
        "sel_gtfs": CheckboxGroup(labels=[os.path.basename(p) for p in gtf_paths], active=active_gtf_indices),
        "spn_height": Spinner(title="Plot Height", low=100, high=2000, step=50, value=config["coverage_and_junctions_profile"]["plot_height"]),
        "chk_xaxis_main": CheckboxGroup(labels=["Show X-Axis"], active=[0] if config["coverage_and_junctions_profile"].get("show_xaxis", False) else []),
        "chk_legend_main": CheckboxGroup(labels=["Show Legend"], active=[0] if config["coverage_and_junctions_profile"].get("show_legend", False) else []),
        "sld_mq": Slider(title="Min Mapping Quality", start=0, end=60, step=1, value=config["coverage_and_junctions_profile"]["min_mapping_quality"]),
        "chk_show_bg": CheckboxGroup(labels=["Show Coverage for MAPQ=0"], active=[0] if config["coverage_and_junctions_profile"].get("coverage_components", {}).get("show_background", False) else []),
        "chk_full_cov": CheckboxGroup(labels=["Show Coverage"], active=[0] if config["coverage_and_junctions_profile"].get("coverage_components", {}).get("show_foreground", False) else []),
        "chk_normalize": CheckboxGroup(labels=["Normalize by Samples"], active=[0] if config["coverage_and_junctions_profile"].get("normalize_by_samples", False) else []),
        "chk_normalize_junctions": CheckboxGroup(labels=["Normalize Junctions"], active=[0] if config["coverage_and_junctions_profile"].get("normalize_junctions_by_samples", True) else []),
        "chk_global_y_range_normalize": CheckboxGroup(labels=["Adapt Global Y-Range to Normalization"], active=[0] if config["coverage_and_junctions_profile"].get("adapt_global_y_range_to_normalization", True) else []),
        "chk_cliffs": CheckboxGroup(labels=["Show Coverage Cliffs"], active=[0]),
        "chk_ambiguity": CheckboxGroup(labels=["Show Ambiguity"], active=[]),
        "sld_ambiguity_factor": Slider(title="Ambiguity Factor", start=1.0, end=10.0, step=0.1, value=config["coverage_and_junctions_profile"]["high_ambiguity_highlighting"]["ambiguity_highlight"].get("threshold", 2.5)),
        "sld_min_reads": Slider(title="Min Reads for Junction", start=1, end=100, step=1, value=1),
        "chk_filter_flanks": CheckboxGroup(labels=["Filter Mismatch/Insertion Flanks"], active=[]),
        "shared_rules_container": column(sizing_mode="stretch_width"),
        "chk_types": CheckboxGroup(labels=["Marked Junctions", "Curated", "Predicted", "Novel"], active=[0, 1, 2, 3]),

        "spn_ann_height": Spinner(title="Annotation Height", low=50, high=1000, step=50, value=config["genome_annotations"]["plot_height"]),
        "chk_xaxis_ann": CheckboxGroup(labels=["Show X-Axis"], active=[0] if config["genome_annotations"].get("show_xaxis", False) else []),
        "chk_markers": CheckboxGroup(labels=["Show Intron Direction"], active=[0] if config["genome_annotations"].get("show_intron_direction", False) else []),
        "chk_cds": CheckboxGroup(labels=["Show CDS Blocks"], active=[0]),
        "chk_ann_filter": CheckboxGroup(labels=["Full Range Annotations"], active=[0] if config["genome_annotations"].get("show_full_range", False) else []),        "chk_unsupported_introns": CheckboxGroup(labels=["Mark Unsupported Introns"], active=[]),
        "spn_max_reads": Spinner(title="Max Extracted Reads", low=10, high=1000, step=10, value=config["full_read_layout"].get("max_extracted_reads", 100)),
        "spn_reads_height": Spinner(title="Reads Plot Height", low=50, high=2000, step=50, value=config["full_read_layout"]["plot_height"]),
        "chk_xaxis_reads": CheckboxGroup(labels=["Show X-Axis"], active=[0] if config["full_read_layout"].get("show_xaxis", False) else []),
        "chk_legend_reads": CheckboxGroup(labels=["Show Legend"], active=[0] if config["full_read_layout"].get("show_legend", False) else []),
        "chk_show_gap_size": CheckboxGroup(labels=["Show Gap Size"], active=[]),
        "chk_show_deletions": CheckboxGroup(labels=["Show Deletion Markers"], active=[0] if config["full_read_layout"].get("show_deletion_markers", True) else []),

        "sel_backend": Select(title="Output Backend", options=["canvas", "svg"], value=config["general"].get("output_backend", "canvas")),
        "chk_fixed_width": CheckboxGroup(labels=["Fixed Plot Width"], active=[0] if config["general"].get("fixed_width", False) else []),
        "spn_plot_width": Spinner(title="Fixed Width (px)", low=400, high=4000, step=100, value=config["general"].get("plot_width", 1200)),

        "chk_crosshair": CheckboxGroup(labels=["Show Sync Crosshair"], active=[0]),
        "chk_crosshair_ext": CheckboxGroup(labels=["Extend Sync Crosshair (Coverage Only)"], active=[]),

        "last_max_y_ann": None,
        "last_ann_state": None,
        "last_data_state": None,
        "last_show_full_cov": False,
        "last_show_bg": False,
        "last_creator_state": None,
        "last_ui_samples": [],
        "last_selected_junctions": [],
        "force_show_cohort_selection": False,

        "div_creator_info": Div(text="", visible=False),
        "btn_add_j": Button(label="Add Selected Junctions to Creator", button_type="success", sizing_mode="stretch_width"),
        "mul_curr_j": MultiSelect(title="Current Junctions (0)", options=[], value=[], size=5, sizing_mode="stretch_width"),
        "btn_remove_j": Button(label="Remove Selected", button_type="warning", sizing_mode="stretch_width"),
        "btn_clear_j": Button(label="Clear All", button_type="danger", sizing_mode="stretch_width"),
        "txt_transcript_id": TextInput(title="New Transcript ID", value="", sizing_mode="stretch_width"),
        "sel_strand": Select(title="Strand", options=["+", "-"], value="+", sizing_mode="stretch_width"),
        "num_t_start": NumericInput(title="Manual Start Position", value=None, sizing_mode="stretch_width"),
        "num_t_end": NumericInput(title="Manual End Position", value=None, sizing_mode="stretch_width"),
        "btn_fetch_start": Button(label="Use View Start", sizing_mode="stretch_width"),
        "btn_fetch_end": Button(label="Use View End", sizing_mode="stretch_width"),
        "btn_zoom_to_creator": Button(label="Zoom to Creator Bounds", button_type="primary", sizing_mode="stretch_width"),
        "btn_create_t": Button(label="Create and Inject Transcript", button_type="success", sizing_mode="stretch_width"),
        "div_download": Div(text="", visible=False),
        "btn_export_json": Button(label="Export Session (JSON)", button_type="primary", sizing_mode="stretch_width"),
        "btn_export_gtf": Button(label="Export as GTF", button_type="primary", sizing_mode="stretch_width"),
        "btn_export_gff3": Button(label="Export as GFF3", button_type="primary", sizing_mode="stretch_width"),
        "file_import_json": FileInput(accept=".json", sizing_mode="stretch_width"),
        "chk_xaxis_creator": CheckboxGroup(labels=["Show X-Axis"], active=[0] if config.get("transcript_creator", {}).get("show_xaxis", False) else []),
        "chk_markers_creator": CheckboxGroup(labels=["Show Intron Direction"], active=[0] if config.get("transcript_creator", {}).get("show_intron_direction", False) else []),
        "chk_autofocus_creator": CheckboxGroup(labels=["Auto-Focus on Selection"], active=[0] if config.get("transcript_creator", {}).get("auto_focus", True) else []),
        "div_transcript_summary": Div(text="""
            <div style="background-color: #f8fbff; border: 1px solid #d0e7ff; border-radius: 6px; padding: 12px; margin-top: 10px; font-family: inherit; color: #2c3e50; line-height: 1.4;">
                <b style="color: #339af0;">Transcript Summary</b><br>
                <span style="color: #868e96; font-style: italic; font-size: 0.9em;">No junctions added.</span>
            </div>
        """),
        "debounce_mq": None,
        "debounce_sample_sel": None,
        "debounce_amb_mq": None,
        "debounce_squeeze": None,
        "selection_updating": False,
        "is_redrawing": False,
        "is_squeezing": False,
        "is_fetching_reads": False,
        "is_sticky_message": False,
        "last_v_step": 1.0,
        "last_total_h": 0.0,
        "user_junctions": [],
        "cohort_selection_trigger": NumericInput(value=0, visible=False),
        "shared_rules_cache": {},
        "src_shared_rules": ColumnDataSource(data=dict(sample=[], curated=[], predicted=[], novel=[])),
        "ds_reads_busy": ColumnDataSource(data=dict(is_busy=[0]))
    }
    return L

def setup_ui(app):
    L = app.L
    config = app.state.get("config", {})
    spacing_val = config.get("general", {}).get("sidebar_plot_spacing", 15)
    # Section Header Labels
    h_data = Div(text="<div style='border-bottom: 2px solid #444; margin-top: 5px; font-weight: bold; color: #333;'>Data Selection</div>", sizing_mode="stretch_width")
    h_coverage = Div(text="<div style='border-bottom: 2px solid #444; margin-top: 15px; font-weight: bold; color: #333;'>Coverage & Junctions Profile</div>", sizing_mode="stretch_width")
    h_cov_visuals = Div(text="<div style='border-bottom: 1px solid #999; margin-top: 10px; font-weight: bold; color: #666;'>└─ Coverage Components</div>", sizing_mode="stretch_width")
    h_junctions = Div(text="<div style='border-bottom: 1px solid #999; margin-top: 10px; font-weight: bold; color: #666;'>└─ Junctions & Splicing</div>", sizing_mode="stretch_width")
    h_ambiguity = Div(text="<div style='border-bottom: 1px solid #999; margin-top: 10px; font-weight: bold; color: #666;'>└─ High Ambiguity Highlighting</div>", sizing_mode="stretch_width")
    h_creator = Div(text="", sizing_mode="stretch_width")
    h_visual_creator = Div(text="<div style='border-bottom: 1px solid #999; margin-top: 10px; font-weight: bold; color: #666;'>└─ Visual Control</div>", sizing_mode="stretch_width")
    h_annotations = Div(text="<div style='border-bottom: 2px solid #444; margin-top: 15px; font-weight: bold; color: #333;'>Genome Annotations</div>", sizing_mode="stretch_width")
    h_reads = Div(text="<div style='border-bottom: 2px solid #444; margin-top: 15px; font-weight: bold; color: #333;'>Full Read Layout</div>", sizing_mode="stretch_width")
    h_general = Div(text="<div style='border-bottom: 2px solid #444; margin-top: 15px; font-weight: bold; color: #333;'>General</div>", sizing_mode="stretch_width")

    tab1_content = column(
        h_data,
        L["sel_gene"], 
        L["btn_show_cohort_selection"],
        L["sel_samples"],
        h_coverage,
        L["spn_height"], 
        L["chk_xaxis_main"], 
        L["chk_legend_main"],
        L["sld_mq"],
        column(h_cov_visuals, L["chk_show_bg"], L["chk_full_cov"], L["chk_normalize"], L["chk_normalize_junctions"], L["chk_global_y_range_normalize"], L["chk_cliffs"], sizing_mode="stretch_width", styles={"margin-left": "15px"}),
        column(h_ambiguity, L["chk_ambiguity"], L["sld_ambiguity_factor"], sizing_mode="stretch_width", styles={"margin-left": "15px"}),
        column(h_junctions, L["sld_min_reads"], L["chk_filter_flanks"], Div(text="<b>Comparative Highlighting Rules:</b>", styles={"margin-top": "5px"}), L["shared_rules_container"], L["chk_types"], sizing_mode="stretch_width", styles={"margin-left": "15px"}),
        h_annotations,
        L["div_gtf_header"],
        L["sel_gtfs"],
        L["spn_ann_height"], 
        L["chk_xaxis_ann"],
        L["chk_markers"], L["chk_cds"], L["chk_ann_filter"], L["chk_unsupported_introns"],
        h_reads,
        L["spn_max_reads"], 
        L["spn_reads_height"], 
        L["chk_xaxis_reads"], 
        L["chk_legend_reads"], 
        L["chk_show_gap_size"], 
        L["chk_show_deletions"],
        h_general,
        L["sel_backend"],
        L["chk_fixed_width"], 
        L["spn_plot_width"],
        L["chk_crosshair"],
        L["chk_crosshair_ext"],
        sizing_mode="stretch_width"
    )

    tab2_content = column(
        h_creator, L["div_creator_info"],
        L["btn_add_j"], L["mul_curr_j"], L["btn_remove_j"], L["btn_clear_j"],
        L["txt_transcript_id"], 
        L["sel_strand"],
        column(L["num_t_start"], L["btn_fetch_start"], sizing_mode="stretch_width", styles={"border": "1px solid #ccc", "padding": "5px", "margin": "1px"}),
        column(L["num_t_end"], L["btn_fetch_end"], sizing_mode="stretch_width", styles={"border": "1px solid #ccc", "padding": "5px", "margin": "1px"}),
        L["btn_zoom_to_creator"], L["btn_create_t"],
        
        Div(text="<div style='border-bottom: 1px solid #999; margin-top: 10px; font-weight: bold; color: #666;'>└─ Export / Import</div>", sizing_mode="stretch_width"),
        L["btn_export_json"],
        L["btn_export_gtf"], 
        L["btn_export_gff3"],
        column(
            Div(text="<b>Re-import Session:</b>", styles={"margin-top": "8px"}),
            L["file_import_json"],
            sizing_mode="stretch_width"
        ),
        L["div_download"],

        column(h_visual_creator, L["chk_xaxis_creator"], L["chk_markers_creator"], L["chk_autofocus_creator"], sizing_mode="stretch_width"),
        L["div_transcript_summary"],
        sizing_mode="stretch_width"
    )

    tabs_list = [
        TabPanel(child=tab1_content, title="Analysis Controls"), 
        TabPanel(child=tab2_content, title="Transcript Editor")
    ]

    L["tabs"] = Tabs(tabs=tabs_list, sizing_mode="stretch_width")
    sidebar = column(L["div_title"], L["div_progress"], L["tabs"], sizing_mode="fixed", width=330, styles={'padding-right': '8px', 'overflow-x': 'hidden'})
    L["plot_column"] = column(L["div_manual"], L["div_config_warnings"], sizing_mode="stretch_width")
    app.layout = row(sidebar, L["plot_column"], sizing_mode="stretch_width", spacing=spacing_val)
    app.doc.title = "lncRNA Seeker Hub"
