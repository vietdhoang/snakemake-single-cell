import fire
import os
import pandas as pd
import scanpy as sc

from anndata import AnnData
from qnorm import quantile_normalize
from scipy.stats import zscore
from typing import Union, List, Tuple

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class Normalize:

    def __init__(self, adata: AnnData, norm_method: str, apply_log: bool = True, 
                 **kwargs) -> None:
        '''Constructor for Normalize class'''
        
        self.adata = adata

        # Convert norm_method (string) into a callable function
        try:
            self.norm_method = getattr(Normalize, norm_method)        
        
        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")
        
        self.norm_method_kwargs = kwargs
        self.apply_log = apply_log


    def run(self) -> None:
        self.norm_method(**self.norm_method_kwargs)

        if self.apply_log:
            self.log()
       
    
    def nonorm(self) -> None:
        '''Don't filter the data. Essentially rename the input file so that
        it still follows the workflow naming convention

        Args:
            path_in: Path to input h5ad file
            path_out: Path where the output h5ad file will be written.
        '''
        pass

    
    def log(self) -> None:
        '''Logarithmize the data'''
        sc.pp.log1p(self.adata)


    def sum(self) -> None:
        '''Normalize by dividing each observation by the sum of their row vector
        '''
        df = self.adata.to_df()
        df.apply(lambda row: row / row.sum())

        self.adata.X = df


    def max(self) -> None:
        '''Normalize by dividing each observation by the maximum value of 
        their row vector
        '''
        df = self.adata.to_df()
        df.apply(lambda row: row / row.max(), axis=1)

        self.adata.X = df


    def median(self) -> None:
        '''Normalize by dividing each observation by the median of their row vector
        '''
        df = self.adata.to_df()
        df.apply(lambda row: row / row.median())

        self.adata.X = df


    def upper_quartile(self) -> None:
        '''Normalize by dividing each observation by the upper quartile 
        of their row vector
        '''
        df = self.adata.to_df()
        df.apply(lambda row: row / row.quantile(0.75))

        self.adata.X = df

    def quantile(self) -> None:
        '''Quantile normalization on each observation (each row will be
        standardized)'''
        self.adata.X = quantile_normalize(self.adata.to_df(), axis=0)


def main(snakemake) -> None:
    # Read in input file 
    adata = sc.read_h5ad(snakemake.input[0])
   
    # Filter the data using the provided filter method
    adata_filtered = Normalize(adata, snakemake.params.norm_method,
                               log=snakemake.params.log).run().adata

    # Write out the file               
    adata_filtered.write(snakemake.output[0])


if __name__ == "__main__":
    if 'snakemake' in globals():
        main(snakemake)
    else:
        fire.Fire()