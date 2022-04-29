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
