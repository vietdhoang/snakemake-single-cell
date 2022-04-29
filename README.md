<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/40422305810c90e2b1b0cbaafc2413bc4530983c/img/dag_one_sample.svg" width=50% height=50%>
</p>

# Snakemake Combinatorial Pipeline

This Snakemake pipeline aims to streamline the process of analyzing single-cell data.   1 
by executing each step in a combinatorial manner. More specifically, if multiple
methods are provided for each stage in the pipeline, this pipeline will automatically 
execute every possible combination  of methods and their parameters. This feature allows 
for effortless exploration of data in a reproducible and scalable manner.

The pipeline performs the following analysis steps:
* Filtering
* Normalization
* Dimensionality Reduction
* Clustering
* Differential Expression Analysis
* Plotting

Below is are some examples of the combinatorial nature of the pipeline. 

Fig. 1 demonstrates how specifying three different normalization methods and three 
different normalization methods will result in 9 different results, each one being one 
possible normalization-dimensionality reduction combination.

<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/40422305810c90e2b1b0cbaafc2413bc4530983c/img/combinatorial_method.svg" width=50% height=50%>
</p>

Fig. 2 demonstrates the possibility of providing multiple parameters to a method. The
pipeline will automatically run every possible combination of parameters for that method.

<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/51ed4eca04327e48a097c6fae04e2d80d9e31bb2/img/combinatorial_param.svg" width=50% height=50%>
</p>

Text
