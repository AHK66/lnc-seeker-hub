# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
import os
import sys
import threading
import logging
import lnc_seeker

# Add current directory to path so we can import our package
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from bokeh.plotting import curdoc
from bokeh.models import Div
from lnc_seeker_bokeh.ui_manager import VisualizerApp
from lnc_seeker_bokeh.pipeline import run_analysis_thread
from lnc_seeker_bokeh.state import verify_environment, load_base_data

doc = curdoc()
# Initialize the app. It will create its own session state.
app = VisualizerApp(doc)
doc.app_instance = app

# Now use the app's state for initialization
state = app.state

# Suppress bokeh.server.protocol.receiver warnings if requested in config
if state.get("config") and state["config"].get("general", {}).get("suppress_bokeh_warnings"):
    logging.getLogger('bokeh.server.protocol.receiver').setLevel(logging.ERROR)

verify_environment(state, lnc_seeker)

if state.get("config") is not None:
    # Kick off backend analysis in background
    threading.Thread(target=run_analysis_thread, args=(state,), daemon=True).start()
else:
    doc.add_root(Div(text="<h1>Failed to initialize environment</h1><p>Check if config.json exists. Check console logs.</p>"))
