# CURRENTLY UNDER DEVELOPMENT
# TODO: Implement the plotting functions for both DimReduce and Cluster

import fire
import os
import scanpy as sc

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
        self.func(self.X, *self.func_args, **self.func_kwargs)
    

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



def pca(path: Union[str, bytes, os.PathLike], 
        path_out: Union[str, bytes, os.PathLike],
        get_plots: bool = True) -> AnnData:
    '''Perform PCA on the data

    TODO:
        Complete documentation by adding more comments towards the end.
        Produce plots that adhere to lab standards.

    Args:
        adata: AnnData object containing the data that PCA will be run on.
    
    Returns:
        AnnData object after PCA
    '''
    # Read in file
    adata = sc.read_h5ad(path)

    # Perform the PCA
    pca = DimReduce(sc.tl.pca, adata, svd_solver='arpack')
    pca.run()
    
    # Access the PCA'ed data
    adata = pca.X

    # Plot PCA results if the option is selected
    if get_plots:
        sc.pl.pca(adata, save=True)
        sc.pl.pca_variance_ratio(adata, log=False, save=True)
    
    # Write out PCA results if the option is selected
    if path_out:
        adata.write(path_out)
    
    
def umap(path: Union[str, bytes, os.PathLike], 
         path_out: Union[str, bytes, os.PathLike],
         get_plots: bool = True) -> AnnData:
    '''Perform upmap on the data and plot results 

    TODO:
        Complete documentation by adding more comments towards the end.
        Write a function that automatically determines the number of principal
        components to use.

    Args:
        adata: AnnData object containing the data for computing the graph
        get_plots: Whether or not to generate graphs associated with this
            procedure.
        out_path: Output path for results file
    
    Returns:
        AnnData object
    '''    

    # Read in file
    adata = sc.read_h5ad(path)

    if not 'X_pca' in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    if not 'neighbors' in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Apply UMAP to the data
    umap = DimReduce(sc.tl.umap, adata)
    umap.run()

    adata = umap.X

    # Plot results if the option is selected.
    if get_plots:
        sc.pl.umap(adata, save=True)

    # Write out UMAP results if the option is selected
    if path_out:
        adata.write(path_out)
    

def tsne(path: Union[str, bytes, os.PathLike], 
         path_out: Union[str, bytes, os.PathLike],
         get_plots: bool = True) -> AnnData:
    '''Perform tsne on the data and plot results 


    Args:
        adata: AnnData object containing the data for computing the graph
        get_plots: Whether or not to generate graphs associated with this
            procedure.
        out_path: Output path for results file
    
    Returns:
        AnnData object
    '''

    # Read in file
    adata = sc.read_h5ad(path)

    if not 'X_pca' in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighbourhood graph. 
    # NOTE: The number of neighbours and components are from the default 
    # settings. This will be amended in the future. 
    if not 'neighbors' in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Apply tSNE to the data
    tsne = DimReduce(sc.tl.tsne, adata)
    tsne.run()

    adata = tsne.X

    # Plot results if the option is selected.
    if get_plots:
        sc.pl.tsne(adata, save=True)

    # Write out tSNE results if the option is selected
    if path_out:
        adata.write(path_out)


def leiden(path: Union[str, bytes, os.PathLike], 
           path_out: Union[str, bytes, os.PathLike],
           get_plots: bool = True) -> AnnData:
    
    # Read in file
    adata = sc.read_h5ad(path)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    leiden = Cluster(sc.tl.leiden, adata)
    leiden.run()

    adata = leiden.X

    if get_plots:
        sc.pl.umap(adata, color='leiden', save=True)
    
    adata.write(path_out)


if __name__ == '__main__':
    fire.Fire()