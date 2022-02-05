import fire
import os
import scanpy as sc
import sys

from typing import Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from scripts.custom.custom_cluster import *

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 1


def leiden(path_in: Union[str, bytes, os.PathLike], 
           path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Cluster the data using the Leiden algorithm

    TODO:
        Ensure all parameters for sc.tl.leiden are determined correctly. 

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
    sc.tl.leiden(adata)
    
    # Write out results
    adata.write(path_out)


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


if __name__ == '__main__':
    fire.Fire()