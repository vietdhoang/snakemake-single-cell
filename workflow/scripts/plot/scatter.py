import altair as alt
import fire
import numpy as np
import os
import pandas as pd
import scanpy as sc
import sys

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_hex
from os.path import dirname
from typing import Final, List, Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, dirname(dirname(dirname(__file__))))
import scripts.common.altair_themes

# Import altair theme from altair_themes.py
alt.themes.register("publish_theme", scripts.common.altair_themes.publish_theme)
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


def get_colour_scheme(colours: List[str], num_colours: int) -> List[str]:
    
    if len(colours) >= num_colours:
        return colours
    
    else:
        cmap = LinearSegmentedColormap.from_list("cmap", colours)
        scheme = cmap(np.linspace(0, 1, num_colours))   

        return [to_hex(c, keep_alpha=False) for c in scheme]


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
            scheme = get_colour_scheme([*SET1.values()], df[label].nunique())

            if label in [*snakemake.config['cluster']]:
                sort = [*map(
                    lambda x: str(x), 
                    sorted(map(lambda x: int(x), df[label].unique()))
                )]
            else:
                sort = sorted(df[label].unique())
            
            chart = alt.Chart(df).mark_circle().encode(
                x=f'{use_rep} 1',
                y=f'{use_rep} 2',
                color=alt.Color(
                    label, 
                    scale=alt.Scale(range=scheme),
                    sort=sort,
                    legend=alt.Legend(title=label.capitalize(), columns=len(sort)//20+1)                    
                ),
                
                tooltip=df.columns.tolist(),
                opacity=alt.condition(selection, alt.value(1), alt.value(0.2))
            ).add_selection(
                selection
            ).properties(
                width=300,
                height=300
            ).interactive()
            
            if label in [*snakemake.config['cluster']]:
                chart.save(f"{path_dir_out}/{prefix}.html")
            else:
                chart.save(
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

           
# if __name__ == "__main__":
#     if 'snakemake' in globals():
        
#         if snakemake.rule == "plot_scatter_labels":
#             path_dir_out = (
#                 f"{snakemake.config['output_dir']}/"
#                 f"{snakemake.wildcards.sample}/figures/labels/"
#             )
#             labels = snakemake.config['labels']        
        
#         elif snakemake.rule == "plot_scatter_dim_reduce":
#             path_dir_out = (
#                 f"{snakemake.config['output_dir']}/"
#                 f"{snakemake.wildcards.sample}/figures/no_labels/"
#             )
#             labels = []
        
#         else:
#             path_dir_out = (
#                 f"{snakemake.config['output_dir']}/"
#                 f"{snakemake.wildcards.sample}/figures/cluster/"
#             )
#             labels = [snakemake.wildcards.c_method]
        
#         prefix = (
#             f"scatter_{snakemake.wildcards.filter_method}_"
#             f"{snakemake.wildcards.norm_method}_"
#             f"{snakemake.wildcards.dr_method}"
#         )
#         use_rep = snakemake.wildcards.dr_method

#         scatter(
#             snakemake.input[0],
#             path_dir_out,
#             prefix,
#             labels,
#             use_rep 
#         )

#     else:
#         fire.Fire()


if __name__ == "__main__":
    if 'snakemake' in globals():
        wildcards = snakemake.wildcards
        if snakemake.rule == "plot_scatter_labels":
            path_dir_out = (
                f"{snakemake.config['output_dir']}/"
                f"{wildcards.sample}/figures/dim_reduce/"
                f"dim_reduce_{wildcards.dr_method}_{wildcards.dr_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
            )
            labels = snakemake.config['labels']        
        
        elif snakemake.rule == "plot_scatter_dim_reduce":
            path_dir_out = (
                f"{snakemake.config['output_dir']}/"
                f"{wildcards.sample}/figures/dim_reduce/"
                f"dim_reduce_{wildcards.dr_method}_{wildcards.dr_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
            )
            labels = []
        
        else:
            path_dir_out = (
                f"{snakemake.config['output_dir']}/"
                f"{wildcards.sample}/figures/cluster/"
                f"cluster_{wildcards.c_method}_{wildcards.c_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                f"dim_reduce_{wildcards.dr_method}_{wildcards.dr_params}/"
            )
            labels = [snakemake.wildcards.c_method]
        
        prefix = "scatter"
        use_rep = snakemake.wildcards.dr_method

        scatter(
            snakemake.input[0],
            path_dir_out,
            prefix,
            labels,
            use_rep 
        )

    else:
        fire.Fire()