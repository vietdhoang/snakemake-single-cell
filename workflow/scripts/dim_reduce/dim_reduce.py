import fire
import scanpy as sc

from anndata import AnnData

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class DimReduce:
    def __init__(self, adata: AnnData, method: str, **kwargs) -> None:
        """Constructor for Filter class"""

        # Convert method (string) into a callable function
        try:
            self.method = getattr(DimReduce, method)

        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")

        self.method_kwargs = kwargs
        self.adata = adata

    def run(self) -> None:
        self.adata = self.method(self.adata, **self.method_kwargs)

    @staticmethod
    def pca(adata: AnnData, n_comps=100, **kwargs) -> AnnData:
        """Perform PCA on the data

        TODO:
            Ensure all parameters for sc.tl.pca are determined correctly.

        Args:
            path_in: AnnData object containing the data that PCA will be run on.
            path_out: Output path for results file
        """

        # Apply tSNE to the data
        sc.tl.pca(adata, n_comps=n_comps, **kwargs)

        return adata

    @staticmethod
    def tsne(adata: AnnData, perplexity=1000, learning_rate=1000, **kwargs) -> AnnData:
        """Perform tsne on the data and plot results

        TODO:
            Write a function that automatically determines the number of principal
            components to use.

        Args:
            path_in: AnnData object containing the data for computing the graph
            path_out: Output path for results file
        """

        # Perform initial PCA to reduce the the computational costs of running
        # t-SNE.
        adata_copy = adata.copy()
        sc.pp.pca(adata_copy, n_comps=50)

        # Apply tSNE to the data
        sc.tl.tsne(
            adata_copy, perplexity=perplexity, learning_rate=learning_rate, **kwargs
        )

        adata.obsm["X_tsne"] = adata_copy.obsm["X_tsne"]

        return adata

    @staticmethod
    def umap(
        adata: AnnData,
        n_neighbours: int = 15,
        min_dist: float = 0.5,
        n_comps: int = 2,
        metric: str = "euclidean",
        **kwargs
    ) -> AnnData:
        """Perform umap on the data and plot results

        TODO:
            Write a function that automatically determines the number of principal
            components to use.

        Args:
            adata: AnnData object containing the data for computing the graph
            path_out: Output path for results file
        """

        # Create neighbourhood graph
        sc.pp.neighbors(
            adata, n_neighbors=n_neighbours, metric=metric, key_added="neighbors_umap"
        )

        # Apply UMAP to the data
        sc.tl.umap(
            adata,
            neighbors_key="neighbors_umap",
            min_dist=min_dist,
            n_components=n_comps,
            **kwargs
        )

        return adata


def main(snakemake) -> None:
    # Read in input file
    adata = sc.read_h5ad(snakemake.input[0])

    # Perform dimensionality reduction with the given method
    adata_dim_reduce = DimReduce(
        adata, snakemake.wildcards.dr_method, **snakemake.params.params
    )
    adata_dim_reduce.run()

    # Write out the file
    adata_dim_reduce.adata.write(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake)
    else:
        fire.Fire()
