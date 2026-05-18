# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
from bokeh.models import Div, Select
from bokeh.layouts import row
from lnc_seeker_bokeh.constants import set_progress_redrawing

class RulesManager:
    def __init__(self, app):
        self.app = app
        self.L = app.L

    def update_rules_ui(self):
        selected_names = self.L["sel_samples"].value
        new_samples, new_curated, new_predicted, new_novel = [], [], [], []
        
        for name in selected_names:
            rules = self.L["shared_rules_cache"].get(name, {"curated": "o (Ignored)", "predicted": "o (Ignored)", "novel": "+ (Present)"})
            new_samples.append(name)
            new_curated.append(rules.get("curated", "o (Ignored)"))
            new_predicted.append(rules.get("predicted", "o (Ignored)"))
            new_novel.append(rules.get("novel", "+ (Present)"))
            self.L["shared_rules_cache"][name] = rules
            
        self.L["src_shared_rules"].data = dict(sample=new_samples, curated=new_curated, predicted=new_predicted, novel=new_novel)

        current_ui_samples = self.L.get("last_ui_samples", [])
        if list(selected_names) != current_ui_samples:
            self.L["last_ui_samples"] = list(selected_names)
            widgets = []
            header = row(
                Div(text="<div style='font-weight: bold; font-size: 0.7em; color: #555;'>Sample</div>", width=120),
                Div(text="<div style='font-weight: bold; font-size: 0.7em; color: #555; transform: rotate(-45deg); transform-origin: bottom left; white-space: nowrap; margin-left: 15px;'>novel</div>", width=45),
                Div(text="<div style='font-weight: bold; font-size: 0.7em; color: #555; transform: rotate(-45deg); transform-origin: bottom left; white-space: nowrap; margin-left: 15px;'>predicted</div>", width=45),
                Div(text="<div style='font-weight: bold; font-size: 0.7em; color: #555; transform: rotate(-45deg); transform-origin: bottom left; white-space: nowrap; margin-left: 15px;'>curated</div>", width=45),
                sizing_mode="stretch_width", styles={"background-color": "#f2f2f2", "border-bottom": "1px solid #ccc", "padding": "2px 4px", "height": "55px", "align-items": "flex-end"}
            )
            widgets.append(header)
            
            for i, name in enumerate(selected_names):
                rules = self.L["shared_rules_cache"].get(name)
                
                def make_cb(s_name, r_type):
                    def cb(attr_cb, old_cb, new_val):
                        self.L["shared_rules_cache"][s_name][r_type] = new_val
                        self.L["is_redrawing"] = True
                        set_progress_redrawing(self.L["div_progress"])
                        # We use add_next_tick_callback to ensure it runs on the Bokeh document thread
                        self.app.doc.add_next_tick_callback(lambda: self.app.on_sample_selection_change(None, None, None))
                    return cb
                
                sel_cur = Select(options=["+ (Present)", "- (Absent)", "o (Ignored)"], value=rules["curated"], width=45, height=26, sizing_mode="fixed")
                sel_cur.on_change('value', make_cb(name, "curated"))
                
                sel_pre = Select(options=["+ (Present)", "- (Absent)", "o (Ignored)"], value=rules["predicted"], width=45, height=26, sizing_mode="fixed")
                sel_pre.on_change('value', make_cb(name, "predicted"))
                
                sel_nov = Select(options=["+ (Present)", "- (Absent)", "o (Ignored)"], value=rules["novel"], width=45, height=26, sizing_mode="fixed")
                sel_nov.on_change('value', make_cb(name, "novel"))
                
                bg_color = "#ffffff" if i % 2 == 0 else "#fafafa"
                display_name = (name[:23] + "...") if len(name) > 23 else name
                label = Div(text=f"<div style='margin-top:2px; font-size: 0.75em; overflow:hidden; text-overflow:ellipsis; white-space: nowrap;' title='{name}'>{display_name}</div>", width=120)
                widgets.append(row(label, sel_nov, sel_pre, sel_cur, sizing_mode="stretch_width", styles={"background-color": bg_color, "border-bottom": "1px solid #eee", "padding": "1px 4px", "align-items": "center"}))
            
            self.L["shared_rules_container"].children = widgets
