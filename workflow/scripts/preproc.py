import altair as alt
import fire
import os
import pandas as pd
import scanpy as sc
import sys

from anndata._core.anndata import AnnData
from typing import List, Union

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import scripts.altair_themes
from scripts.custom.custom_qc import *

# Import altair theme from altair_themes.py
alt.themes.register("publish_theme", scripts.altair_themes.publish_theme)
alt.themes.enable("publish_theme")

# Set verbosity to 1 which means only print out errors
sc.settings.verbosity = 1


def mtx_to_h5ad(path_in: Union[str, bytes, os.PathLike],
                path_out: Union[str, bytes, os.PathLike],
                path_label: Union[str, bytes, os.PathLike] = None,
                labels: List[str] = None) -> None:
    '''Convert a count matrix in '.mtx' format to '.h5ad' format.

    TODO:
        May have to change printed error message into a log entry for a error
        log
    
    Args:
        path: Path to directory containing barcode.tsv, genes.tsv, and matrix.mtx
        path_out: Path where the '.h5ad' file will be written.
        prefix: Prefix for barcode.tsv.gz, genes.tsv.gz and matrix.mtx.gz files.
            ex. 'subject1_barcode.tsv.gz' -> prefix is 'subject1_'
        path_label: Optionally provide a path to a label file which will be
            incorporated into the AnnData for analysis.
        labels: List of labels from the label file that will be included in
            the analysis
    '''

    adata = sc.read_10x_mtx(path_in, var_names='gene_symbols')
    
    if path_label:
        df_labels = pd.read_csv(path_label, index_col=0)
        
        for col_name in labels:
            adata.obs[col_name] = df_labels[col_name].astype('category')

    if os.path.splitext(path_out)[1] == '.h5ad':          
        adata.write(path_out)
    
    else:
        print("Error: Output file needs to have a '.h5ad' extension",
              file=sys.stderr)              
        exit(1)


def plot_qc_metrics(adata: AnnData) -> None:
    '''Plot QC metrics.

    Args:
        adata: AnnData object containing the data to be plotted.
        path_dir_out: path where output figures will be written.

    '''
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                                       log1p=False, inplace=True)

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


def noqc(path_in: Union[str, bytes, os.PathLike],
         path_out: Union[str, bytes, os.PathLike]) -> None:
    '''Don't perfom QC on the data. Essentially rename the input file so that
    it still follows the workflow naming convention

    Args:
        path_in: Path to input h5ad file
        path_out: Path where the output h5ad file will be written.
    '''

    # Read in input file 
    adata = sc.read_h5ad(path_in)

    # Validate output path and write out the file
    if os.path.splitext(path_out)[1] == '.h5ad':                
        adata.write(path_out)


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


if __name__ == '__main__':
    fire.Fire()