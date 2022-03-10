import fire
import os
import scanpy as sc

from typing import Union

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


def nonorm(path_in: Union[str, bytes, os.PathLike],
           path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Don't perfom normalization on the data. Essentially rename the input file
    so that it still follows the workflow naming convention

    Args:
        path_in: Path to input h5ad file
        path_out: Path where the output h5ad file will be written.
    '''

    # Read in input file 
    adata = sc.read_h5ad(path_in)

    # Validate output path and write out the file
    if os.path.splitext(path_out)[1] == '.h5ad':                
        adata.write(path_out)


if __name__ == "__main__":
    if 'snakemake' in globals():
        nonorm(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()