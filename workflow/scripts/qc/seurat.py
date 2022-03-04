import fire
import os
import scanpy as sc

from typing import Union


def seurat(path_in: Union[str, bytes, os.PathLike],
           path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform QC on the data using the same steps as what is outlined in the
    Seurat tutorial (Satija et al., 2015). The steps are outlined in this paper:
    https://www.nature.com/articles/nbt.3154

    Args:
        path_in: Path to h5ad file containing the data that will be qc'ed
        path_out: Path where the output h5ad file will be written.
    '''
    
    # Read in input file 
    adata = sc.read_h5ad(path_in)

    # plot_qc_metrics(adata)

    # Used scanpy's Seurat recipe for the QC 
    sc.pp.recipe_seurat(adata)
    
    # Validate output ptah and write out the qc'ed file
    if os.path.splitext(path_out)[1] == '.h5ad':                
        adata.write(path_out)


if __name__ == "__main__":
    if 'snakemake' in globals():
        seurat(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()