# CURRENTLY UNDER DEVELOPMENT
# TODO: Implement the plotting functions for both DimReduce and Cluster

import os

from numpy import ndarray
from scipy.sparse import spmatrix 
from anndata._core.anndata import AnnData
from typing import Callable, Union


class DataMatrixTransform:
    
    def __init__(self, func: Callable, X: Union[AnnData, ndarray, spmatrix],
                 *func_args, **func_kwargs) -> None:        
        
        self.func = func
        self.func_args = func_args
        self.func_kwargs = func_kwargs
        self.X = X
    
    def run(self) -> None:
        self.func(self.X, self.func_args, self.func_kwargs)
    

class DimReduce(DataMatrixTransform):

    def __init__(self, func: Callable, X: Union[AnnData, ndarray, spmatrix], 
                 *func_args, **func_kwargs) -> None:
        super().__init__(func, X, *func_args, **func_kwargs)
    
    
    def plot(self, path_dir_out: Union[str, bytes, os.PathLike]) -> None:
        pass


class Cluster(DataMatrixTransform):


    def __init__(self, func: Callable, X: Union[AnnData, ndarray, spmatrix], 
                 *func_args, **func_kwargs) -> None:
        super().__init__(func, X, *func_args, **func_kwargs)

    
    def plot(self, path_dir_out: Union[str, bytes, os.PathLike]) -> None:
        pass