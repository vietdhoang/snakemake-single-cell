import fire
import pandas as pd
import scanpy as sc

from anndata import AnnData
from sklearn.cluster import KMeans

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class Cluster:
    def __init__(self, adata: AnnData, method: str, **kwargs) -> None:
        """Constructor for Cluster class"""

        # Convert method (string) into a callable function
        try:
            self.method = getattr(Cluster, method)

        except AttributeError:
            raise NotImplementedError("{method_name} has not been implemented")

        self.method_kwargs = kwargs
        self.adata = adata

    def run(self) -> None:
        """Run the culestring algorithm"""
        self.adata = self.method(self.adata, **self.method_kwargs)

    @staticmethod
    def leiden(
        adata: AnnData,
        n_neighbours=15,
        use_rep: str = "X",
        metric="euclidean",
        **kwargs,
    ) -> AnnData:
        """Cluster the data using the Leiden algorithm

        Args:
            adata: AnnData object containing the data that will be clustered.
            n_neighbours: Number of neighbours to construct the neighbourhood graph.
            use_rep: Which data space to perform clustering on (i.e X_pca -> cluster
                on PCA space. X -> cluster on using the original data)
            metric: Which distance metric to use to determine similarity

        Returns:
            AnnData containing cluster labels
        """

        # Create neighbourhood graph
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbours,
            metric=metric,
            use_rep=use_rep,
            key_added="neighbors_leiden",
        )

        # Cluster the data
        sc.tl.leiden(adata, neighbors_key="neighbors_leiden", **kwargs)

        return adata

    @staticmethod
    def louvain(
        adata: AnnData,
        n_neighbours: int = 15,
        use_rep: str = "X",
        metric: str = "euclidean",
        **kwargs,
    ) -> AnnData:
        """Cluster the data using the Louvain algorithm

        Args:
            adata: AnnData object containing the data that will be clustered.
            n_neighbours: Number of neighbours to construct the neighbourhood graph.
            use_rep: Which data space to perform clustering on (i.e X_pca -> cluster
                on PCA space. X -> cluster on using the original data)
            metric: Which distance metric to use to determine similarity

        Returns:
            AnnData containing cluster labels
        """

        # Create neighbourhood graph
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbours,
            metric=metric,
            use_rep=use_rep,
            key_added="neighbors_louvain",
        )

        # Cluster the data
        sc.tl.louvain(adata, neighbors_key="neighbors_louvain", **kwargs)

        return adata

    @staticmethod
    def kmeans(
        adata: AnnData, n_clusters: int = 8, use_rep: str = "X", **kwargs
    ) -> AnnData:
        """Cluster the data using the K-Means algorithm. This algorithm uses scikit
        learn's implementation

        Args:
            adata: AnnData object containing the data that will be clustered.
            n_clusters: number of clusters to create from the data
            use_rep: Which data space to perform clustering on (i.e X_pca -> cluster
                on PCA space. X -> cluster on using the original data)

        Returns:
            AnnData containing cluster labels
        """

        if use_rep == "X":
            mat = adata.X.toarray()
        else:
            mat = adata.obsm[use_rep]

        labels = KMeans(n_clusters, **kwargs).fit(mat).labels_
        adata.obs["kmeans"] = pd.Categorical(labels)

        return adata


def main(snakemake) -> None:
    # Read in input file
    adata = sc.read_h5ad(snakemake.input[0])

    # Perform dimensionality reduction with the given method
    adata_dim_reduce = Cluster(
        adata,
        snakemake.wildcards.c_method,
        use_rep=snakemake.params.representation,
        **snakemake.params.params,
    )
    adata_dim_reduce.run()

    # Write out the file
    adata_dim_reduce.adata.write(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake)
    else:
        fire.Fire()
