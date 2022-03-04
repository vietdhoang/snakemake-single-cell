import fire
import os
import scanpy as sc

from typing import Union


def tutorial(path_in: Union[str, bytes, os.PathLike],
             path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform QC the EXACT same way that was done in the scanpy tutorial:
    https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
    This should only be run on the 3K pbmc dataset.

    Args:
        path: Path to h5ad file containing the data that will be qc'ed
        path_out: Path where the output h5ad file will be written.
    '''
    
    sc.settings.verbosity = 0

    # Read in input file 
    adata = sc.read_h5ad(path_in)

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                               log1p=False, inplace=True)
    
    # Remove cells that have too many mitochondrial genes expressed 
    # or too many total counts 
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    # Total-count normalize (library-size correct) the data matrix X to 
    # 10,000 reads per cell, so that counts become comparable among cells.
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logorithmize the data
    sc.pp.log1p(adata)

    # Identify highly variable genes and keep them
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    # Regress out effects of total counts per cell and the percentage of 
    # mitochondrial genes expressed. Scale the data to unit variance.
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # Scale each gene to unit variance. 
    # Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)
    
    # Validate output ptah and write out the qc'ed file
    if os.path.splitext(path_out)[1] == '.h5ad':                
        adata.write(path_out)


if __name__ == "__main__":
    if 'snakemake' in globals():
        tutorial(snakemake.input[0], snakemake.output[0])    
    else:
        fire.Fire()