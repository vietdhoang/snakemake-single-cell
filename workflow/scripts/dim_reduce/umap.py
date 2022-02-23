import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 0


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


if __name__ == "__main__":
    if 'snakemake' in globals():
        umap(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()