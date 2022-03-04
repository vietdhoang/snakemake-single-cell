import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 0


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


if __name__ == "__main__":
    if 'snakemake' in globals():
        neighbourhood(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()