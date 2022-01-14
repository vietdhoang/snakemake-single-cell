import os
import scanpy as sc
from anndata._core.anndata import AnnData
from typing import Union


def plot_qc_metrics(adata: AnnData, 
                    path_dir_out: Union[str, bytes, os.PathLike]) -> None:
    '''Plot QC metrics.

    Args:
        adata: AnnData object containing the data to be plotted.
        path_dir_out: path where output figures will be written.

    '''
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata = sc.calculate_qc_metrics

    sc.pl.violin(
        adata, 
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4, 
        multi_panel=True,
        save='_computed_qc_measures.png'
    )

    # Show mitochondrial gene expression
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', 
                  save='_mit_gene_expr.png')

    # Show total gene count per cell
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  save='_total_gene_count.png')


def qc_seurat(path: Union[str, bytes, os.PathLike],
              path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Perform QC on the data using the same steps as what is outlined in the
    Seurat tutorial (Satija et al., 2015). The steps are outlined in this paper:
    https://www.nature.com/articles/nbt.3154

    Args:
        path: Path to h5ad file containing the data that will be qc'ed
        path_out: Path where the output h5ad file will be written.
    '''
    
    # Read in input file 
    adata = sc.read_h5ad(path)

    # Used scanpy's Seurat recipe for the QC 
    sc.pp.recipe_seurat(adata)
    
    # Validate output ptah and write out the qc'ed file
    if os.path.splitext(path_out)[1] == '.h5ad':                
        adata.write(path_out)