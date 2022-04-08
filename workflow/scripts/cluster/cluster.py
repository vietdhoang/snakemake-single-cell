import fire
import os 
import scanpy as sc

from anndata import AnnData

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class Cluster:

    def __init__(self, adata: AnnData, method: str, **kwargs) -> None:
        '''Constructor for Filter class'''
        
        # Convert method (string) into a callable function
        try:
            self.method = getattr(Cluster, method)        
        
        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")
        
        self.method_kwargs = kwargs
        self.adata = adata


    def run(self) -> None:
        self.adata = self.method(self.adata, **self.method_kwargs)
    

    @staticmethod
    def leiden(adata: AnnData, n_neighbours=15, use_rep: str = 'X',
               metric='euclidean', **kwargs) -> AnnData:
        '''Cluster the data using the Leiden algorithm

        TODO:
            Ensure all parameters for sc.tl.leiden are determined correctly. 

        Args:
            path_in: AnnData object containing the data that will be clustered.
            path_out: Output path for results file
        '''

        # Create neighbourhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbours, metric=metric,
                        use_rep=use_rep, key_added="neighbors_leiden")

        # Cluster the data
        sc.tl.leiden(adata, neighbors_key="neighbors_leiden")
        
        return adata
    

    @staticmethod
    def louvain(adata: AnnData, n_neighbours: int = 15, use_rep: str = 'X',
                metric: str = 'euclidean', **kwargs) -> AnnData:
        '''Cluster the data using the Louvain algorithm

        TODO:
            Ensure all parameters for sc.tl.louvain are determined correctly. 

        Args:
            path_in: AnnData object containing the data that will be clustered.
            path_out: Output path for results file
        '''    
        
        # Create neighbourhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbours, metric=metric,
                        use_rep=use_rep, key_added="neighbors_louvain")

        # Cluster the data
        sc.tl.leiden(adata, neighbors_key="neighbors_louvain")
        
        return adata


def main(snakemake) -> None:
    # Read in input file 
    adata = sc.read_h5ad(snakemake.input[0])

    # Perform dimensionality reduction with the given method
    adata_dim_reduce = Cluster(adata, 
                               snakemake.wildcards.c_method,
                               use_rep=snakemake.params.representation,
                               **snakemake.params.params)
    adata_dim_reduce.run()

    # Write out the file               
    adata_dim_reduce.adata.write(snakemake.output[0])


if __name__ == "__main__":
    if 'snakemake' in globals():
        main(snakemake)
    else:
        fire.Fire()