import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 0


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


if __name__ == "__main__":
    if 'snakemake' in globals():
        tsne(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()