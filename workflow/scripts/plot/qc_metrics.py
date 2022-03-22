import altair as alt
import fire
import os
import pandas as pd
import scanpy as sc
import sys

from os.path import join as pjoin
from os.path import dirname
from typing import Final, List, Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, dirname(dirname(dirname(__file__))))
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

# Define the set1 Vega color scheme for reference
QC_METRIC: Final[dict] = {
    'n_genes_by_counts':    'Number of genes expressed by the cell',
    'n_cells_by_counts':    'Number of cells that express the gene',
    'total_counts':         'Total UMI Counts',
    'pct_counts_mt':        'Percentage of genes that are mitochondrial'
}


def violin_plot_qc(df, path_dir, basename):
    
    for col in df:
        alt.Chart(df).transform_density(
            col,
            as_=[QC_METRIC[col], 'Density'],
        ).mark_area().encode(
            x=alt.X(QC_METRIC[col], type='quantitative'),
            y='Density:Q'
        ).properties(
            width=300,
            height=300
        ).interactive().save(pjoin(path_dir, f"{basename}_{col}.html"))
    


def main(snakemake):

    # Read in input file 
    adata = sc.read_h5ad(snakemake.input[0])

    # Plot QC metrics for cells
    df_cells = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']]
    violin_plot_qc(df_cells, snakemake.params.dir_out, "density_cells")
    
    # Plot QC metrics for genes
    df_cells = adata.obs[['n_cells_by_counts', 'total_counts']]
    violin_plot_qc(df_cells, snakemake.params.dir_out, "density_genes")




if __name__ == "__main__":

    if 'snakemake' in globals():
        main(snakemake)
