# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
# Visual progress styles and HTML helpers
PROGRESS_STYLES = """
<style>
@keyframes pulse-yellow {
  0% { box-shadow: 0 0 0 0px rgba(255, 193, 7, 0.4); }
  70% { box-shadow: 0 0 0 10px rgba(255, 193, 7, 0); }
  100% { box-shadow: 0 0 0 0px rgba(255, 193, 7, 0); }
}
@keyframes pulse-blue {
  0% { box-shadow: 0 0 0 0px rgba(33, 150, 243, 0.4); }
  70% { box-shadow: 0 0 0 10px rgba(33, 150, 243, 0); }
  100% { box-shadow: 0 0 0 0px rgba(33, 150, 243, 0); }
}
@keyframes pulse-red {
  0% { box-shadow: 0 0 0 0px rgba(220, 53, 69, 0.4); }
  70% { box-shadow: 0 0 0 10px rgba(220, 53, 69, 0); }
  100% { box-shadow: 0 0 0 0px rgba(220, 53, 69, 0); }
}
@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}
.ui-spinner {
  border: 3px solid #f3f3f3;
  border-top: 3px solid #e67e22;
  border-radius: 50%;
  width: 16px;
  height: 16px;
  animation: spin 1s linear infinite;
  display: inline-block;
  vertical-align: middle;
  margin-right: 8px;
}
.progress-card {
    padding: 10px;
    border-radius: 8px;
    margin-bottom: 10px;
    box-shadow: 0 2px 5px rgba(0,0,0,0.1);
    font-family: sans-serif;
    height: 60px;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    justify-content: center;
    gap: 4px;
}
</style>
"""

def set_progress_in_progress(div):
    set_progress_message(div, "Analysis in progress...", True)

def set_progress_message(div, message, show_spinner=False):
    spinner_html = '<div class="ui-spinner"></div>' if show_spinner else ""
    new_text = PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #fff3cd; border: 1px solid #ffeeba; animation: pulse-yellow 2s infinite;">
       <div>
           {spinner_html}
           <span style="color: #856404; font-weight: bold; font-size: 13px;">{message}</span>
       </div>
    </div>
    """
    if div.text != new_text:
        div.text = new_text

def set_progress_fail(div, message):
    new_text = PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #f8d7da; border: 1px solid #f5c6cb;">
       <div>
           <span style="color: #721c24; font-weight: bold; font-size: 13px;">❌ {message}</span>
       </div>
    </div>
    """
    if div.text != new_text:
        div.text = new_text

def set_progress_complete(div):
    new_text = PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #d4edda; border: 1px solid #c3e6cb;">
       <div>
           <span style="color: #155724; font-weight: bold; font-size: 13px;">✅ Scan Complete.</span>
       </div>
    </div>
    """
    if div.text != new_text:
        div.text = new_text

def set_progress_success(div, message):
    new_text = PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #d4edda; border: 1px solid #c3e6cb;">
       <div>
           <span style="color: #155724; font-weight: bold; font-size: 13px;">✅ {message}</span>
       </div>
    </div>
    """
    if div.text != new_text:
        div.text = new_text

def clear_progress(div):
    # Sets an empty card to maintain height and prevent UI "jumping"
    new_text = PROGRESS_STYLES + '<div class="progress-card" style="box-shadow: none; border: none; background: transparent;"></div>'
    if div.text != new_text:
        div.text = new_text

def set_progress_redrawing(div):
    new_text = PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #e8f4fd; border: 1px solid #b8daff; animation: pulse-blue 2s infinite;">
       <div>
           <div class="ui-spinner" style="border-top: 3px solid #2196F3;"></div>
           <span style="color: #004085; font-weight: bold; font-size: 13px;">Redrawing Profiles...</span>
       </div>
       <div style="font-size: 11px; color: #004085; margin-left: 24px;">Applying Comparative Highlighting Rules</div>
    </div>
    """
    if div.text != new_text:
        div.text = new_text

# Define custom Red-Gray-Blue diverging palette for change visualization
RedGrayBlue11 = [
    '#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', # Distinct Blues (Negative)
    '#444444',                                           # Neutral (Dark Gray)
    '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d'  # Distinct Reds (Positive)
]

def get_progress_html(stage_text, percent, prog_str, bar_color, text_color):
    return PROGRESS_STYLES + f"""
    <div class="progress-card" style="background-color: #f8f9fa; border: 1px solid #dee2e6;">
        <div style="font-weight: bold; color: #333; font-size: 13px;">
            <div class="ui-spinner"></div>{stage_text}...
        </div>
        <div style="width: 100%; background-color: #e9ecef; border-radius: 10px; height: 18px; position: relative; overflow: hidden; border: 1px solid #ced4da;">
            <div style="width: {percent}%; height: 100%; background-color: {bar_color}; transition: width 0.4s ease-out;"></div>
            <div style="position: absolute; top: 0; width: 100%; text-align: center; font-size: 11px; line-height: 18px; font-weight: bold; color: {text_color};">
                {prog_str}
            </div>
        </div>
    </div>
    """
