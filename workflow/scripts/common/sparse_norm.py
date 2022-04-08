import numpy as np
import scipy

from anndata import AnnData
from functools import wraps
from numpy.typing import ArrayLike, NDArray
from qnorm import quantile_normalize
from typing import Callable, NewType


def _check_anndata_sparsity(func: Callable) -> Callable:
    '''Decorator for sprase normalization functions. Checks if the data in 
    the X attribute of an AnnData is a sparse matrix. If X is not sparse,
    then convert it into a sparse matrix
    '''

    @wraps(func)
    def wrapper(adata: AnnData, **kwargs) -> AnnData:
        if not scipy.sparse.issparse(adata.X):
            adata.X = scipy.sparse.csr_array(adata.X)
        
        return func(adata, **kwargs)

    return wrapper
    
    
def _apply_rows(func: Callable, adata: AnnData, **kwargs) -> AnnData:
    '''Apply a function row-wise to the X attribute of an AnnData object'''
    
    mat = adata.X
    
    for i in range(mat.shape[0]):
        row = mat.getrow(i)
        row_transformed_data = func(row.data, **kwargs)
        mat.data[mat.indptr[i]: mat.indptr[i+1]] = row_transformed_data
    
    return adata
    

@_check_anndata_sparsity
def log1p_sparse(adata: AnnData) -> AnnData:
    '''Apply natural logarithm to the data'''
    adata.X = adata.X.log1p()
    return adata
    

@_check_anndata_sparsity
def sum_norm_sparse(adata: AnnData) -> AnnData:
    '''Normalize by dividing by each row by the row-sum'''
        
    def sum_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.sum(arr)

    adata = _apply_rows(sum_norm_1D, adata)

    return adata
    

@_check_anndata_sparsity
def max_norm_sparse(adata: AnnData) -> AnnData:
    '''Normalize by dividing by each row by the row-max'''

    def max_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.max(arr)

    adata = _apply_rows(max_norm_1D, adata)

    return adata


@_check_anndata_sparsity
def percentile_norm_sparse(adata: AnnData, 
                           percentile: float = 50) -> AnnData:
    '''Normalize by dividing by each row by the row-percentile'''
        
    def percentile_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.percentile(arr, percentile, interpolation='midpoint')
        
    adata = _apply_rows(percentile_norm_1D, adata)

    return adata
    

@_check_anndata_sparsity
def quantile_norm_sparse(adata: AnnData) -> AnnData:
    '''Quantile normalization on each observation (each row will be
    standardized)'''

    mat = np.array(adata.X.todense())
    adata.X = scipy.sparse.csr_matrix(quantile_normalize(mat, axis=0))
    
    return adata

