import altair as alt
import fire
import os
import pandas as pd
import scanpy as sc
import sys

from typing import Final, List, Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import scripts.altair_themes

# Import altair theme from altair_themes.py
alt.themes.register("publish_theme", scripts.altair_themes.publish_theme)
alt.themes.enable("publish_theme")

# Define the set1 Vega color scheme for reference
SET1: Final[dict] = {
    'red':      '#e41a1c',
    'blue':     '#377eb8',
    'green':    '#4daf4a',
    'purple':   '#984ea3',
    'orange':   '#ff7f00',
    'yellow':   '#ffff33',
    'brown':    '#a65628',
    'pink':     '#f781bf',
    'gray':     '#999999'
}


def scatter(path_in: Union[str, bytes, os.PathLike], 
            path_dir_out: Union[str, bytes, os.PathLike],
            prefix: str, 
            labels: List[str],
            use_rep: str) -> None:
    
    # Read in input file 
    adata = sc.read_h5ad(path_in)
    
    df = pd.DataFrame({
        'index': adata.obs_names.astype('string'),
        'observation': adata.obs_names.astype('string'),
        f'{use_rep} 1': adata.obsm[f'X_{use_rep}'][:, 0],
        f'{use_rep} 2': adata.obsm[f'X_{use_rep}'][:, 1]
    }).set_index('index')

    # Add in all the observation annotations to display while hovering over the
    # points in the interactive plot
    for annot in adata.obs_keys():
        df[annot] = adata.obs[annot]

    if labels:
        for label in labels:
            selection = alt.selection_multi(fields=[label], bind='legend')
            alt.Chart(df).mark_circle().encode(
                x=f'{use_rep} 1',
                y=f'{use_rep} 2',
                color=alt.Color(label, scale=alt.Scale(scheme='set1')),
                tooltip=df.columns.tolist(),
                opacity=alt.condition(selection, alt.value(1), alt.value(0.2))
            ).add_selection(
                selection
            ).properties(
                width=300,
                height=300
            ).interactive().save(
                f"{path_dir_out}/{prefix}_{label}.html"
            )
    
    else:
        alt.Chart(df).mark_circle(opacity=1).encode(
            x=f'{use_rep} 1',
            y=f'{use_rep} 2',
            color=alt.value(SET1['blue']),
            tooltip=df.columns.tolist()
        ).properties(
            width=300,
            height=300
        ).interactive().save(
            f"{path_dir_out}/{prefix}.html"
        )

           
if __name__ == '__main__':
    fire.Fire()