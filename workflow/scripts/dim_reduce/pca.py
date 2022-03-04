import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 0


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


if __name__ == "__main__":
    if 'snakemake' in globals():
        pca(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()