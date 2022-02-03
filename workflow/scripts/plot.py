import altair as alt
import fire
import os
import pandas as pd
import scanpy as sc
import sys

from typing import Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import scripts.altair_themes
from scripts.custom.custom_dim_reduce import *

# Import altair theme from altair_themes.py
alt.themes.register("publish_theme", scripts.altair_themes.publish_theme)
alt.themes.enable("publish_theme")
