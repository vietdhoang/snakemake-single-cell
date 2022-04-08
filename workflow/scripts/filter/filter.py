import fire
import pandas as pd
import scanpy as sc
import sys

from anndata import AnnData
from os.path import dirname
from typing import Callable, Tuple

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, dirname(dirname(dirname(__file__))))
import scripts.common.filter_cutoffs as filter_cutoffs

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class Filter:

    def __init__(self, adata: AnnData, method: str, **kwargs) -> None:
        '''Constructor for Filter class'''
        
        # Convert method (string) into a callable function
        try:
            self.method = getattr(Filter, method)        
        
        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")
        
        self.method_kwargs = kwargs
        self.adata = adata


    def run(self) -> None:
        self.adata = Filter.basic_filter(self.adata)
        self.adata = self.method(self.adata, **self.method_kwargs)


    @staticmethod
    def _query(df: pd.DataFrame, 
               col_name: str, 
               lower_cutoff: float = None, 
               upper_cutoff: float = None) -> list[str]:
        
        if not upper_cutoff and not lower_cutoff:
            raise ValueError("At least one of upper_bound or lower_bound must be True")        
        
        elif upper_cutoff and lower_cutoff:            
            condition = f"{col_name} < {lower_cutoff} or {col_name} > {upper_cutoff}"        
        
        elif lower_cutoff:
            condition = f"{col_name} < {lower_cutoff}"
        
        else:
            condition = f"{col_name} > {upper_cutoff}"
        
        return df.query(condition).index.tolist()        
    
    
    @staticmethod
    def _blacklist(get_cutoffs: Callable, adata: AnnData) -> Tuple[list[float]]:
        # Observation and variable blacklist (ones that will be filtered out)
        obs_bl = set()
        var_bl = set()

        # Total cell UMI count filtering (observation filtering)
        # Filter out cells that have too high or too low total UMI counts across
        # all genes
        lo, hi = get_cutoffs(adata.obs['total_counts'])
        obs_bl = obs_bl.union(set(Filter._query(adata.obs, 
                                                'total_counts',
                                                lower_cutoff=lo, 
                                                upper_cutoff=hi)))
        
        # Number of genes filtering (observation filtering)
        # Filter out cells that have too few genes expressed.
        lo, hi = get_cutoffs(adata.obs['n_genes_by_counts'])
        obs_bl = obs_bl.union(set(Filter._query(adata.obs, 
                                                'n_genes_by_counts',
                                                lower_cutoff=lo,
                                                upper_cutoff=hi)))
        
        # Percentage mitochondrial gene filtering (observation filtering)
        # Filter out cells that have too high or too low percentage of total
        # counts that are mitochondrial
        lo, hi = get_cutoffs(adata.obs['pct_counts_mt'])
        obs_bl = obs_bl.union(set(Filter._query(adata.obs, 
                                                'pct_counts_mt',
                                                lower_cutoff=lo,
                                                upper_cutoff=hi)))
        
        # Total gene UMI count filtering (variable filtering)
        # Filter out genes that have too low total UMI counts across all cells
        var_bl = var_bl.union(set(Filter._query(adata.var, 
                                                'total_counts',
                                                lower_cutoff=lo)))
        
        # Number of cells filtering (variable filtering)
        # Filter out genes that have too few cells that express them
        lo, hi = get_cutoffs(adata.var['n_cells_by_counts'])
        var_bl = var_bl.union(set(Filter._query(adata.var, 
                                                'n_cells_by_counts',
                                                lower_cutoff=lo)))
        return (list(obs_bl), list(var_bl))


    @staticmethod
    def nofilter(adata: AnnData) -> AnnData:
        '''Don't filter the data. Essentially rename the input file so that
        it still follows the workflow naming convention

        Args:
            path_in: Path to input h5ad file
            path_out: Path where the output h5ad file will be written.
        '''
        return adata


    @staticmethod
    def basic_filter(adata: AnnData, 
                     min_genes: int = 1, 
                     min_cells: int = 1) -> AnnData:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
    
        return adata


    @staticmethod
    def IQR(adata: AnnData, k: float = 1.5) -> AnnData:
        '''Interquartile range outlier filtering.'''
        
        def get_cutoffs(adata: AnnData) -> Tuple[float]:
            return filter_cutoffs.get_iqr_cutoffs(adata, k=k)
        
        obs_blacklist, var_blacklist = Filter._blacklist(get_cutoffs, adata)
        
        adata = adata[~adata.obs_names.isin(obs_blacklist), 
                      ~adata.var_names.isin(var_blacklist)]
        
        return adata


    @staticmethod
    def MAD(adata: AnnData, threshold: float = 2.5) -> AnnData:
        '''Mean absolute deviation filtering'''
        
        def get_cutoffs(adata: AnnData) -> Tuple[float]:
            return filter_cutoffs.get_mad_cutoffs(adata, 
                                                  threshold=threshold)
        
        obs_blacklist, var_blacklist = Filter._blacklist(get_cutoffs, adata)
        
        adata = adata[~adata.obs_names.isin(obs_blacklist), 
                      ~adata.var_names.isin(var_blacklist)]
        
        return adata


    @staticmethod
    def zscore(adata: AnnData, threshold: float = 1.96) -> AnnData:
        '''z-score filtering'''
        
        def get_cutoffs(adata: AnnData) -> Tuple[float]:
            return filter_cutoffs.get_zscore_cutoffs(adata, 
                                                     threshold=threshold)
        
        obs_blacklist, var_blacklist = Filter._blacklist(get_cutoffs, adata)
        
        adata = adata[~adata.obs_names.isin(obs_blacklist), 
                      ~adata.var_names.isin(var_blacklist)]
        
        return adata


def main(snakemake) -> None:
    # Read in input file 
    adata = sc.read_h5ad(snakemake.input[0])

    # Filter the data using the provided filter method
    adata_filtered = Filter(adata, 
                            snakemake.wildcards.filter_method,
                            **snakemake.params.params)
    adata_filtered.run()

    # Write out the file               
    adata_filtered.adata.write(snakemake.output[0])


if __name__ == "__main__":
    if 'snakemake' in globals():
        main(snakemake)
    else:
        fire.Fire()