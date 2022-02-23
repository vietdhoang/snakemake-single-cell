import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


def louvain(path_in: Union[str, bytes, os.PathLike], 
            path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Cluster the data using the Louvain algorithm

    TODO:
        Ensure all parameters for sc.tl.louvain are determined correctly. 

    Args:
        path_in: AnnData object containing the data that will be clustered.
        path_out: Output path for results file
    '''    
    
    # Read in file
    adata = sc.read_h5ad(path_in)

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    # Cluster the data
    sc.tl.louvain(adata)
    
    # Write out results
    adata.write(path_out)


if __name__ == "__main__":
    if 'snakemake' in globals():
        louvain(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()