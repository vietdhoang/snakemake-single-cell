import fire
import scanpy as sc
import sys

from anndata import AnnData
from os.path import dirname

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, dirname(dirname(dirname(__file__))))
import scripts.common.sparse_norm as sparse_norm


class Norm:
    def __init__(
        self, adata: AnnData, method: str, apply_log: bool = True, **kwargs
    ) -> None:
        """Constructor for Normalize class"""

        self.adata = adata

        # Convert norm_method (string) into a callable function
        try:
            self.method = getattr(Norm, method)

        except AttributeError:
            raise NotImplementedError(f"{method} has not been implemented")

        self.method_kwargs = kwargs
        self.apply_log = False if self.method == "log1p" else apply_log

    def run(self) -> None:
        self.adata = self.method(self.adata, **self.method_kwargs)

        if self.apply_log:
            self.adata = Norm.log1p(self.adata)

    @staticmethod
    def nonorm(adata: AnnData) -> AnnData:
        """Don't filter the data. Essentially rename the input file so that
        it still follows the workflow naming convention

        Args:
            path_in: Path to input h5ad file
            path_out: Path where the output h5ad file will be written.
        """
        return adata

    @staticmethod
    def log1p(adata: AnnData) -> AnnData:
        """Logarithmize the data"""
        return sparse_norm.log1p_sparse(adata)

    @staticmethod
    def sum(adata: AnnData) -> AnnData:
        """Normalize by dividing each observation by the sum of their row vector"""
        return sparse_norm.sum_norm_sparse(adata)

    @staticmethod
    def max(adata: AnnData) -> AnnData:
        """Normalize by dividing each observation by the maximum value of
        their row vector
        """
        return sparse_norm.max_norm_sparse(adata)

    @staticmethod
    def median(adata: AnnData) -> AnnData:
        """Normalize by dividing each observation by the median of their
        row vector
        """
        return sparse_norm.percentile_norm_sparse(adata, percentile=50)

    @staticmethod
    def upper_quartile(adata: AnnData) -> AnnData:
        """Normalize by dividing each observation by the upper quartile
        of their row vector
        """
        return sparse_norm.percentile_norm_sparse(adata, percentile=75)

    @staticmethod
    def quantile(adata: AnnData) -> AnnData:
        """Quantile normalization on each observation (each row will be
        standardized)"""
        return sparse_norm.quantile_norm_sparse(adata)


def main(snakemake) -> None:
    # Read in input file
    adata = sc.read_h5ad(snakemake.input[0])

    # Filter the data using the provided filter method
    norm = Norm(adata, snakemake.wildcards.norm_method, **snakemake.params.params)
    norm.run()

    adata_norm = norm.adata

    # Write out the file
    adata_norm.write(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake)
    else:
        fire.Fire()
