import fire
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.stats

from anndata import AnnData
from numpy.typing import ArrayLike
from statsmodels.stats.multitest import multipletests
from typing import List, Tuple

# Set verbosity to 0 which means only print out errors
sc.settings.verbosity = 0


class StatTest:
    @staticmethod
    def test(test: str, x1: ArrayLike, x2: ArrayLike, **kwargs) -> Tuple[float]:
        # Convert method (string) into a callable function
        try:
            if test == "stat_test":
                raise ValueError("'stat_test' cannot be used as a test.")

            func = getattr(StatTest, test)

        except AttributeError:
            raise NotImplementedError("{test} has not been implemented")

        return func(x1, x2, **kwargs)

    @staticmethod
    def mannwhitneyu(x1: ArrayLike, x2: ArrayLike, **kwargs) -> Tuple[float]:
        score, pval = scipy.stats.mannwhitneyu(x1, x2, **kwargs)
        return (score, pval)

    @staticmethod
    def ttest(x1: ArrayLike, x2: ArrayLike, **kwargs) -> Tuple[float]:
        score, pval = scipy.stats.ttest_ind(x1, x2, **kwargs)
        return (score, pval)


def multi_column_isin(df, d):
    # get rid of the index to avoid index alignment
    booleans = [df[col].isin(arr).array for col, arr in d.items()]
    return np.logical_and.reduce(booleans)


def diff_exp_group1_vs_group2(
    adata: AnnData,
    group1: dict,
    group2: dict,
    drop_zeros: bool = False,
    use_raw: bool = False,
    alpha: float = 0.05,
    stat_test: str = "mannwhitneyu",
    multiple_comparisons="fdr_bh",
) -> pd.DataFrame:

    if use_raw:
        adata = adata.raw

    adata1 = adata[multi_column_isin(adata.obs, group1)]
    adata2 = adata[multi_column_isin(adata.obs, group2)]

    dict_results = {
        "name": [],
        "log2foldchange": [],
        "score": [],
        "pvals": [],
        "pvals_adj": [],
    }

    for i in range(adata.shape[1]):

        if drop_zeros:
            col1 = adata1[:, i].X.toarray().flatten()
            col2 = adata2[:, i].X.toarray().flatten()
        else:
            col1 = adata1[:, i].X.data
            col2 = adata2[:, i].X.data

        if len(col1) == 0 or len(col2) == 0:
            continue

        score, p = StatTest.test(stat_test, col1, col2)

        dict_results["name"].append(adata.var_names[i])
        dict_results["log2foldchange"].append(np.log2(np.mean(col1) / np.mean(col2)))
        dict_results["score"].append(score)
        dict_results["pvals"].append(p)

    func_multi_comparisons = lambda x: multipletests(x, alpha, multiple_comparisons)[1]
    dict_results["pvals_adj"] = func_multi_comparisons(dict_results["pvals"])

    df_results = pd.DataFrame(dict_results).set_index("name")
    df_results.sort_values(by=["log2foldchange"], ascending=False, inplace=True)

    return df_results


def main(snakemake) -> None:
    # Read in input file
    adata = sc.read_h5ad(snakemake.input[0])
    params = snakemake.params.params

    if "group1" in params and "group2" in params:
        adata.uns["diff_exp_group1_vs_group2"] = diff_exp_group1_vs_group2(
            adata,
            params["group1"],
            params["group2"],
            **{k: v for k, v in params.items() if k not in ["group1", "group2"]},
        )
    elif "one_vs_rest" in params:
        category_name = params["one_vs_rest"]
        categories = set(adata.obs[category_name].unique())
        for one in categories:
            rest = list(categories - {one})

            group1 = (
                {category_name: [one], **params["filters"]}
                if params["filters"]
                else {category_name: [one]}
            )

            group2 = (
                {category_name: rest, **params["filters"]}
                if params["filters"]
                else {category_name: rest}
            )

            adata.uns[f"diff_exp_{one}_vs_rest"] = diff_exp_group1_vs_group2(
                adata,
                group1,
                group2,
                **{
                    k: v
                    for k, v in params.items()
                    if k not in ["one_vs_rest", "filters"]
                },
            )

    # Write out the file
    adata.write(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake)
    else:
        fire.Fire()
