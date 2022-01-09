import os
import pandas as pd
from pandas.core.algorithms import quantile
import scanpy as sc
import sys

from anndata._core.anndata import AnnData
from typing import Union, NoReturn


def qc(adata: AnnData, get_plots: bool = True) -> AnnData:
    '''Perform quality control

    TODO:
        Complete documentation by adding more comments towards the end.

    Args:
        adata: AnnData object containing the data to be qc'd
    
    Returns:
        AnnData object that has been qc'd
    '''    

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                               log1p=False, inplace=True)

    # Show qc plots if get_plots is set to True (default is False)
    if get_plots:
        
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
    
    # Perform the filtering via slicing.
    # Remove cells that have too many mitochondrial genes expressed or have too
    # many total counts.
    # 1.5x interquartile range
    quartile1 = adata.obs.n_genes_by_counts.quantile(0.25)
    quartile3 = adata.obs.n_genes_by_counts.quantile(0.75)
    
    adata = adata[adata.obs.n_genes_by_counts > 1.5*quartile1, :]
    adata = adata[adata.obs.n_genes_by_counts < 1.5*quartile3, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125,
                                max_mean=3, min_disp=0.5)
    
    if get_plots:
        sc.pl.highly_variable_genes(adata, save='.png')
    
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    return adata