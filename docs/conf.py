from __future__ import annotations

from datetime import datetime

project = "prideQC"
author = "PRIDE Team and contributors"
copyright = f"{datetime.now().year}, PRIDE Team and contributors"
release = "0.2.0"

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    # Maintainer-only figure generators and scientific captions are not public pages.
    "figure_source/**",
    "Thumbs.db",
    ".DS_Store",
    # Repository-local development / benchmarking notes are intentionally
    # excluded from the public end-user documentation site.
    "codon-pride-accessions.md",
    "codon-pride-file-array.md",
    "development-quality.md",
    "hpc.md",
    "mass-error-historical-comparison.md",
    "mass-shift-scout.md",
    "migration.md",
    "slurm-benchmark.md",
    "validation.md",
]

html_theme = "furo"
html_title = "prideQC documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/PRIDE-Archive/prideQC/",
    "source_branch": "main",
    "source_directory": "docs/",
    "footer_icons": [],
    "top_of_page_buttons": ["view"],
    "announcement": (
        "prideQC computes quality-control metrics from raw mass-spectrometry data "
        "and uses that evidence to refine SDRF metadata."
    ),
}

html_logo = None
html_favicon = None
html_show_sourcelink = False
html_permalinks_icon = "#"
myst_heading_anchors = 3
