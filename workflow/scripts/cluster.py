import os
import scanpy as sc
from anndata._core.anndata import AnnData
from typing import Union

def leiden(path: Union[str, bytes, os.PathLike], 
           path_out: Union[str, bytes, os.PathLike],
           get_plots: bool = True) -> AnnData:
    
    # Read in file
    adata = sc.read_h5ad(path)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    sc.tl.leiden(adata)

    if get_plots:
        sc.pl.umap(adata, color='leiden', save=True)
    
    adata.write(path_out)