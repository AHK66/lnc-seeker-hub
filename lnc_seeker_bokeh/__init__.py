# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
"""Lightweight lnc_seeker_bokeh package initializer.

Avoid importing submodules at package import time because several
submodules depend on the compiled `lnc_seeker` extension which may not
be available during development or in CI before the extension is built.

Import submodules explicitly (e.g. ``from lnc_seeker_bokeh import state`` or
``import lnc_seeker_bokeh.state``) where needed so imports remain predictable.
"""

__all__ = [
	"state",
	"data_utils",
	"pipeline",
	"constants",
	"shared_data",
]
