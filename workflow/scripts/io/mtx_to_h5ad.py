import fire
import os
import pandas as pd
import scanpy as sc
import sys

from typing import List, Union


def _get_qc_metrics(adata) -> None:
    # Add an annotation labeling mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Calculate QC metrics. These will be used for filtering
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    return adata


def mtx_to_h5ad(
    path_in: Union[str, bytes, os.PathLike],
    path_out: Union[str, bytes, os.PathLike],
    path_label: Union[str, bytes, os.PathLike] = None,
    labels: List[str] = None,
) -> None:
    """Convert a count matrix in '.mtx' format to '.h5ad' format.

    TODO:
        May have to change printed error message into a log entry for a error
        log

    Args:
        path: Path to directory containing barcode.tsv, genes.tsv, and matrix.mtx
        path_out: Path where the '.h5ad' file will be written.
        prefix: Prefix for barcode.tsv.gz, genes.tsv.gz and matrix.mtx.gz files.
            ex. 'subject1_barcode.tsv.gz' -> prefix is 'subject1_'
        path_label: Optionally provide a path to a label file which will be
            incorporated into the AnnData for analysis.
        labels: List of labels from the label file that will be included in
            the analysis
    """

    adata = sc.read_10x_mtx(path_in, var_names="gene_symbols")

    if path_label:
        df_labels = pd.read_csv(path_label, index_col=0)

        for col_name in labels:
            adata.obs[col_name] = df_labels[col_name].astype("category")

    adata = _get_qc_metrics(adata)

    if os.path.splitext(path_out)[1] == ".h5ad":
        adata.write(path_out)

    else:
        print("Error: Output file needs to have a '.h5ad' extension", file=sys.stderr)
        exit(1)


if __name__ == "__main__":

    if (
        "snakemake" in globals()
        and "label_path" in snakemake.config["inputs"][snakemake.wildcards.sample]
        and snakemake.config["inputs"][snakemake.wildcards.sample]
    ):

        mtx_to_h5ad(
            snakemake.input[0],
            snakemake.output[0],
            path_label=snakemake.config["inputs"][snakemake.wildcards.sample][
                "label_path"
            ],
            labels=snakemake.config["labels"],
        )

    elif "snakemake" in globals():
        mtx_to_h5ad(snakemake.input[0], snakemake.output[0])

    else:
        fire.Fire()
