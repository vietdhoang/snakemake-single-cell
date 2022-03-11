import fire
import os
import pandas as pd
import scanpy as sc

from anndata import AnnData
from scipy.stats import zscore
from typing import Union, List, Tuple

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class Filter:

    def __init__(self, adata: AnnData, filter_method: str, **kwargs) -> None:
        '''Constructor for Filter class'''
        
        # Convert filter_method (string) into a callable function
        try:
            self.filter_method = getattr(Filter, filter_method)        
        
        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")
        
        self.filter_method_kwargs = kwargs

        # Add an annotation labeling mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')

        # Calculate QC metrics. These will be used for filtering
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                                   log1p=False, inplace=True)
        
        self.adata = adata


    def run(self) -> None:
        self.filter_method(**self.filter_method_kwargs)
    
    
    @staticmethod
    def _create_blacklist(df: pd.DataFrame, 
                          col_name: str, 
                          lower_cutoff: float = None, 
                          upper_cutoff: float = None) -> List[str]:
        
        if not upper_cutoff and not lower_cutoff:
            raise ValueError("At least one of upper_bound or lower_bound must be True")        
        elif upper_cutoff and lower_cutoff:            
            condition = f"{col_name} < {lower_cutoff} or {col_name} > {upper_cutoff}"        
        elif lower_cutoff:
            condition = f"{col_name} < {lower_cutoff}"
        else:
            condition = f"{col_name} > {upper_cutoff}"
        
        return df.query(condition).index.tolist()        
    
    
    def nofilter(self) -> None:
        '''Don't filter the data. Essentially rename the input file so that
        it still follows the workflow naming convention

        Args:
            path_in: Path to input h5ad file
            path_out: Path where the output h5ad file will be written.
        '''
        pass

    
    def IQR(self, k: float = 1.5) -> None:
        '''Interquartile range outlier filtering.'''

        # Observation and variable blacklist (ones that will be filtered out)
        obs_bl = {}
        var_bl = {}

        # Total cell UMI count filtering (observation filtering)
        # Filter out cells that have too high or too low total UMI counts across
        # all genes
        q1, q3 = self.adata.obs['total_counts'].quantile([0.25, 0.75])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'total_counts',
                                                            lower_cutoff=k*q1,
                                                            upper_cutoff=k*q3)))
        
        # Number of genes filtering (observation filtering)
        # Filter out cells that have too few genes expressed.
        q1, q3 = self.adata.obs['total_genes_by_counts'].quantile([0.25, 0.75])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'total_genes_by_counts',
                                                            lower_cutoff=k*q1)))

        # Percentage mitochondrial gene filtering (observation filtering)
        # Filter out cells that have too high or too low percentage of total
        # counts that are mitochondrial
        q1, q3 = self.adata.obs['pct_counts_mt'].quantile([0.25, 0.75])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'pct_counts_mt',
                                                            lower_cutoff=k*q1,
                                                            upper_cutoff=k*q3)))
        
        # Total gene UMI count filtering (variable filtering)
        # Filter out genes that have too low total UMI counts across all cells
        q1, q3 = self.adata.obs['total_counts'].quantile([0.25, 0.75])
        var_bl = var_bl.union(set(Filter._create_blacklist(self.adata.var, 
                                                           'total_counts',
                                                            lower_cutoff=k*q1)))
        
        # Number of cells filtering (variable filtering)
        # Filter out genes that have too few cells that express them
        q1, q3 = self.adata.obs['n_cells_by_counts'].quantile([0.25, 0.75])
        var_bl = var_bl.union(set(Filter._create_blacklist(self.adata.var, 
                                                           'n_cells_by_counts',
                                                            lower_cutoff=k*q1)))
        
        self.adata = self.adata[~self.adata.obs_names.isin(list(obs_bl)), 
                                ~self.adata.var_names.isin(list(var_bl))]


    def MAD(self, m: float = 2.5) -> None:
        '''Mean absolute deviation filtering'''
        pass


    def zscore(self, z: float = 1.96) -> None:
        '''z-score filtering'''
        
        def get_cutoffs(s: pd.Series) -> Tuple(float):
            std = s.std()
            mean = s.mean()
            lower_cutoff = -z*std + mean
            upper_cutoff = -z*std + mean

            return (lower_cutoff, upper_cutoff)

        
        # Observation and variable blacklist (ones that will be filtered out)
        obs_bl = {}
        var_bl = {}

        # Total cell UMI count filtering (observation filtering)
        # Filter out cells that have too high or too low total UMI counts across
        # all genes
        lo, hi = get_cutoffs(self.adata.obs['total_counts'])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'total_counts',
                                                            lower_cutoff=lo,
                                                            upper_cutoff=hi)))
        
        # Number of genes filtering (observation filtering)
        # Filter out cells that have too few genes expressed.
        lo, hi = get_cutoffs(self.adata.obs['total_genes_by_counts'])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'total_genes_by_counts',
                                                            lower_cutoff=lo,
                                                            upper_cutoff=hi)))

        # Percentage mitochondrial gene filtering (observation filtering)
        # Filter out cells that have too high or too low percentage of total
        # counts that are mitochondrial
        lo, hi = get_cutoffs(self.adata.obs['pct_counts_mt'])
        obs_bl = obs_bl.union(set(Filter._create_blacklist(self.adata.obs, 
                                                           'pct_counts_mt',
                                                            lower_cutoff=lo,
                                                            upper_cutoff=hi)))
        
        # Total gene UMI count filtering (variable filtering)
        # Filter out genes that have too low total UMI counts across all cells
        lo, hi = get_cutoffs(self.adata.var['total_counts'])
        var_bl = var_bl.union(set(Filter._create_blacklist(self.adata.var, 
                                                           'total_counts',
                                                            lower_cutoff=lo)))
        
        # Number of cells filtering (variable filtering)
        # Filter out genes that have too few cells that express them
        lo, hi = get_cutoffs(self.adata.var['n_cells_by_counts'])
        var_bl = var_bl.union(set(Filter._create_blacklist(self.adata.var, 
                                                           'n_cells_by_counts',
                                                            lower_cutoff=lo)))
        
        self.adata = self.adata[~self.adata.obs_names.isin(list(obs_bl)), 
                                ~self.adata.var_names.isin(list(var_bl))]


def main(snakemake) -> None:
    # Read in input file 
    adata = sc.read_h5ad(snakemake.input[0])
   
    # Filter the data using the provided filter method
    adata_filtered = Filter(adata, snakemake.params.filter_method).run().adata

    # Write out the file               
    adata_filtered.write(snakemake.output[0])


if __name__ == "__main__":
    if 'snakemake' in globals():
        main(snakemake)
    else:
        fire.Fire()