import os
import scanpy as sc
import sys

from anndata._core.anndata import AnnData
from typing import Union


def read_count_matrix(path: Union[str, bytes, os.PathLike], 
                      prefix: str = None) -> AnnData:
    '''Read in count matrix

    TODO:
        May have to change printed error message into a log entry for a error
        log

    Args:
        path: Path to directory cotaining barcode.tsv, genes.tsv, and matrix.mtx
            OR path to .h5ad file. 
    
    Returns:
        AnnData object
    '''

    if os.path.isdir(path):
        adata = sc.read_10x_mtx(path, var_names='gene_symbols', prefix=prefix)

    elif os.path.splitext(path)[1] == '.h5ad':
        adata = sc.read_10x_h5(path)
    
    else:
        print("Error: Unrecognized count matrix input", file=sys.stderr)
        exit(1)
    
    return adata


def mtx_to_h5ad(path: Union[str, bytes, os.PathLike],
                path_out: Union[str, bytes, os.PathLike],
                prefix: str = None) -> None:
    '''Convert a count matrix in '.mtx' format to '.h5ad' format.

    TODO:
        May have to change printed error message into a log entry for a error
        log
    
    Args:
        path: Path to directory containing barcode.tsv, genes.tsv, and matrix.mtx
        path_out: Path where the '.h5ad' file will be written.
        prefix: Prefix for barcode.tsv.gz, genes.tsv.gz and matrix.mtx.gz files.
            ex. 'subject1_barcode.tsv.gz' -> prefix is 'subject1_'
    '''
    
    adata = read_count_matrix(path, prefix=prefix)
    
    if os.path.splitext(path_out)[1] == '.h5ad':          
        adata.write(path_out)
    
    else:
        print("Error: Output file needs to have a '.h5ad' extension",
              file=sys.stderr)              
        exit(1)