import fire
import os
import pandas as pd
import scanpy as sc
import scipy

from typing import Union

def h5ad_to_csv(path_in: Union[str, bytes, os.PathLike],
                path_mtx_out: Union[str, bytes, os.PathLike],
                path_label_out: Union[str, bytes, os.PathLike]) -> None:
    
    adata = sc.read_h5ad(path_in)
    X = adata.X.T.todense() if scipy.sparse.issparse(adata.X) else adata.X.T
    df_mat = pd.DataFrame(X, columns=adata.obs_names, index=adata.var_names)
    df_mat.to_csv(path_mtx_out)

    df_label = pd.DataFrame({'item': adata.obs_names, 'label': adata.obs['label']})
    df_label.to_csv(path_label_out, index=False)


if __name__ == "__main__":

    if 'snakemake' in globals():
        h5ad_to_csv(snakemake.input[0], snakemake.output[0], snakemake.output[1])
    else:
        fire.Fire()