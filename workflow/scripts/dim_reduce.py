import altair as alt
import fire
import os
import pandas as pd
import scanpy as sc
import sys

from typing import Union


def neighbourhood(path_in: Union[str, bytes, os.PathLike], 
                  path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Compute the neighbourhood graph from the data

    TODO:
        Ensure all parameters for sc.pp.neighborhood are determined correctly. 

    Args:
        path_in: AnnData object containing the data.
        path_out: Output path for results file
    ''' 

    # Read in file
    adata = sc.read_h5ad(path_in)

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # Write out PCA results if the option is selected
    adata.write(path_out)# Write out PCA results if the option is selected


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
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    
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
    
    # Apply tSNE to the data
    sc.tl.tsne(adata)
    
    # Write out tSNE results if the option is selected
    adata.write(path_out)


if __name__ == '__main__':
    fire.Fire()