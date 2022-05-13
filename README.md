<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/40422305810c90e2b1b0cbaafc2413bc4530983c/img/dag_one_sample.svg" width=50% height=50%>
</p>

# Snakemake Combinatorial Pipeline

This Snakemake pipeline aims to streamline the process of analyzing single-cell data 
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

Below is are some examples of the combinatorial nature of the pipeline. The dataset
used to generate these figures contain 1500 B-cells, 1500 cytotoxic T-cells, and
1500 helper T-cells. The data were retrieved from 
Zheng et al. *Nature Communications*, 2017.

Fig. 1 demonstrates how specifying three different normalization methods and three 
different normalization methods will result in 9 different results, each one being one 
possible normalization-dimensionality reduction combination.

<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/40422305810c90e2b1b0cbaafc2413bc4530983c/img/combinatorial_method.svg" width=75% height=50%>
</p>

Fig. 2 demonstrates the possibility of providing multiple parameters to a method. The
pipeline will automatically run every possible combination of parameters for that method.

<p align="center">
  <img src="https://github.com/vietdhoang/snakemake-single-cell/blob/51ed4eca04327e48a097c6fae04e2d80d9e31bb2/img/combinatorial_param.svg" width=75% height=50%>
</p>

## Installation and Setup
### Dependencies
* [Conda](https://docs.conda.io/en/latest/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (>=7.4)
* [too-many-cells](https://gregoryschwartz.github.io/too-many-cells/)

### Installation
Clone the repo:
```
git clone https://github.com/vietdhoang/snakemake-single-cell.git
```
Done! Before using the pipeline, make sure to activate the Snakemake environment:
```
conda activate Snakemake
```
### (Optional) Setup for offline use
This post-installation setup step is useful for execution in a HPC environment where internet 
access may not be accessible. The goal is to download all the necessary libraries so that the
pipeline will work even without internet (for example in a compute node).

Activate the Snakemake environment.
```
conda activate snakemake
```
Enter into the cloned repo:
```
cd <PATH/TO/REPO>
```
Run the command below. It create all the necessary conda environments for this
pipeline. The environments will be found in `snakemake-single-cell/.snakemake/conda
```
snakemake create_all_conda_envs -c --use-conda --conda-create-envs-only
```
## Basic Usage
