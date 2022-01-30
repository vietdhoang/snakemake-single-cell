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


def pca(path_in: Union[str, bytes, os.PathLike], 
        path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform PCA on the data

    TODO:
        Ensure all parameters for sc.tl.pca are determined correctly. 

    Args:
        path_in: AnnData object containing the data that PCA will be run on.
        path_out: Output path for results file
    '''
    
    # Read in file
    adata = sc.read_h5ad(path_in)

    # Perform the PCA
    sc.tl.pca(adata, n_comps=100, svd_solver='arpack')
    
    # Write out PCA results if the option is selected
    adata.write(path_out)
    
    
def umap(path_in: Union[str, bytes, os.PathLike], 
         path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform umap on the data and plot results 

    TODO:
        Write a function that automatically determines the number of principal
        components to use.

    Args:
        adata: AnnData object containing the data for computing the graph
        path_out: Output path for results file
    '''    

    # Read in file
    adata = sc.read_h5ad(path_in)

    if not 'X_pca' in adata.obsm:
        sc.tl.pca(adata, n_comps=100, svd_solver='arpack')

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    if not 'neighbors' in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Apply UMAP to the data
    sc.tl.umap(adata)

    # Write out UMAP results
    adata.write(path_out)
    

def tsne(path_in: Union[str, bytes, os.PathLike], 
         path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform tsne on the data and plot results 

    TODO:
        Write a function that automatically determines the number of principal
        components to use.

    Args:
        path_in: AnnData object containing the data for computing the graph
        path_out: Output path for results file
    '''

    # Read in file
    adata = sc.read_h5ad(path_in)

    if not 'X_pca' in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    if not 'neighbors' in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Apply tSNE to the data
    sc.tl.tsne(adata)
    
    # Write out tSNE results if the option is selected
    adata.write(path_out)


if __name__ == '__main__':
    fire.Fire()