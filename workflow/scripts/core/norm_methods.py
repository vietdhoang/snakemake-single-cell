import numpy as np
import scipy

from functools import wraps
from numpy.typing import ArrayLike, NDArray
from qnorm import quantile_normalize
from typing import Callable, NewType

SparseMatrix = NewType("SparseMatrix", scipy.sparse.base.spmatrix)


def _validate_sparse_matrix(func: Callable) -> Callable:
    '''Decorator for sprase normalization functians'''

    @wraps(func)
    def wrapper(mat, **kwargs):
        if scipy.sparse.issparse(mat):
            return func(mat, **kwargs)
        else:
            raise ValueError("Argument must be a Scipy sparse matrix!!")

    return wrapper
    
    
def _apply_sparse(func: Callable, mat: SparseMatrix, **kwargs) -> SparseMatrix:
    for i in range(mat.shape[0]):
        row = mat.getrow(i)
        row_transformed_data = func(row.data, **kwargs)
        mat.data[mat.indptr[i]: mat.indptr[i+1]] = row_transformed_data
    
    return mat
    

@_validate_sparse_matrix
def log1p_sparse(mat: SparseMatrix) -> SparseMatrix:
    '''Apply natural logarithm to the data'''
    return mat.log1p()
    

@_validate_sparse_matrix
def sum_norm_sparse(mat: SparseMatrix) -> SparseMatrix:
    '''Normalize by dividing by each row by the row-sum'''
        
    def sum_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.sum(arr)
        
    return _apply_sparse(sum_norm_1D, mat)
    

@_validate_sparse_matrix
def max_norm_sparse(mat: SparseMatrix) -> SparseMatrix:
    '''Normalize by dividing by each row by the row-max'''

    def max_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.max(arr)

    return _apply_sparse(max_norm_1D, mat)


@_validate_sparse_matrix
def percentile_norm_sparse(mat: SparseMatrix, 
                           percentile: float = 50) -> SparseMatrix:
    '''Normalize by dividing by each row by the row-percentile'''
        
    def percentile_norm_1D(arr: ArrayLike) -> NDArray:
        return arr / np.percentile(arr, percentile, method='midpoint')
        
    return _apply_sparse(percentile_norm_1D, mat)
    

@_validate_sparse_matrix
def quantile_norm_sparse(mat: SparseMatrix) -> SparseMatrix:
    '''Quantile normalization on each observation (each row will be
    standardized)'''
    mat = mat.todense()
    return scipy.sparse.csr_matrix(quantile_normalize(mat, axis=0))

