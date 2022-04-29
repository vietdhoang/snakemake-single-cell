import copy
import itertools
import os
import pandas as pd
import pathlib
import re

from pathlib import Path
from snakemake.utils import Paramspace
from typing import Final, List, Union, Tuple


def get_merge_samples_input(wildcards) -> List[str]:
    """Determines the input files for the rule merge_samples

    Args:
        wildcards: A snakemake object that holds all the wildcards.
    
    Returns:
        A list of required input files
    """
    # Go through the list of samples and add the normalized and filtered version of 
    # that sample to the list of desired inputs
    samples = []
    for sample_name in config["inputs"]:
        samples.append(
            (
                f"{config['output_dir']}/{sample_name}/norm/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                "mtx.h5ad"
            )
        )

    return samples


def get_maketree_input(wildcards):
    """Determines the input files for the rule tmc_maketree

    Args:
        wildcards: A snakemake object that holds all the wildcards.
    
    Returns:
        A list of required input files
    """

    inputs = []

    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.norm_method == "tmc":
        #  If the sample's name is not "merged"
        if wildcards.sample == "merged":
            # Append the path of all inputs
            # The path to a concatonated label file is required as well.
            for input in config["inputs"]:
                inputs.append(config["inputs"][input]["data_path"])

            # Check if the user provided label files by checking
            # if a label file exists for the first input (i.e input_1)
            first_input = [*config["inputs"].keys()][0]
            if exists("label_path", config["inputs"][first_input]):
                inputs.append(
                    (
                        f"{config['output_dir']}/merged/too-many-cells/"
                        f"{wildcards.maketree}/merged_labels.csv"
                    )
                )
        #  If the sample's name is not "merged"
        else:
            inputs.append(config["inputs"][wildcards.sample]["data_path"])

            # Add the label file if provided.
            if exists("label_path", config["inputs"][wildcards.sample]):
                inputs.append(config["inputs"][wildcards.sample]["label_path"])

    # If the user decides to use this pipeline's qc methods instead of
    # too-many-cells' filtering and normalization
    else:
        # Append qc'd matrix csv file to inputs
        inputs.append(
            (
                f"{config['output_dir']}/{wildcards.sample}/norm/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                f"h5ad2csv/mtx.csv"
            )
        )

        # Append corresponding label file for that csv file if it exists
        # If the sample is 'merged' and a label file exists for the first input,
        # then it must exist for all other inputs (user must satisfy this)
        # If the sample is not 'merged' (i.e anything else), the check if the
        # sample has a label file. If it does, add it to the list of inputs.
        first_input = [*config["inputs"].keys()][0]
        if wildcards.sample == "merged":
            if exists("label_path", config["inputs"][first_input]):
                inputs.append(
                    (
                        f"{config['output_dir']}/{wildcards.sample}/norm/"
                        f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                        f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                        f"h5ad2csv/label.csv"
                    )
                )

        elif exists("label_path", config["inputs"][wildcards.sample]):
            inputs.append(
                (
                    f"{config['output_dir']}/{wildcards.sample}/norm/"
                    f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                    f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                    f"h5ad2csv/label.csv"
                )
            )

    return inputs


def get_maketree_params(wildcards):
    """Parses the options for maketree.

    Args:
        wildcards: A snakemake object that holds all the wildcards.
    
    Returns:
        A string representing all the options selected for too-many-cells maketree.
        These options are in the correct format for too-many-cell's command line tool.
        For example, if this function outputs the string s, then too-many-cells will
        be executed as such:

        too-many-cells make-tree s

        As an example, s could be:        
        s = "--matrix-path input --labels-file labels.csv --output out > clusters.csv"
    """

    # Get the list of user-selected options from the config file
    options = config["too-many-cells"][wildcards.maketree]["options"]

    inputs = []

    # Get the first input (input_1)
    first_input = [*config["inputs"].keys()][0]

    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.norm_method == "tmc":
        # If the sample's name is not "merged"
        if wildcards.sample == "merged":
            for input in config["inputs"]:
                inputs.append(f"--matrix-path {config['inputs'][input]['data_path']}")

            # Check if the user provided label files by checking
            # if a label file exists for the first input
            if exists("label_path", config["inputs"][first_input]):
                inputs.append(
                    (
                        f"--labels-file {config['output_dir']}/merged/"
                        f"too-many-cells/{wildcards.maketree}/merged_labels.csv"
                    )
                )

        # If the sample's name is not "merged"
        else:
            # Append qc'd matrix csv file to inputs
            inputs.append(
                f"--matrix-path {config['inputs'][wildcards.sample]['data_path']}"
            )

            # Add the label file if provided.
            if exists("label_path", config["inputs"][wildcards.sample]):
                inputs.append(
                    f"--labels-file {config['inputs'][wildcards.sample]['label_path']}"
                )

    # If the user decides to use this pipeline's qc methods instead of
    # too-many-cells' filtering and normalization
    else:
        # Append qc'd matrix csv file to inputs
        inputs.append(
            (
                f"--matrix-path "
                f"{config['output_dir']}/{wildcards.sample}/norm/"
                f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                f"h5ad2csv/mtx.csv"
            )
        )

        # Append corresponding label file for that csv file (if it exists).
        # If the sample is 'merged' and a label file exists for the first input,
        # then include the --labels-file option for the 'merged' sample
        # If the sample is not 'merged' (i.e anything else), the check if the
        # sample has a label file. If it does, add it to the list of inputs.
        if wildcards.sample == "merged":
            if exists("label_path", config["inputs"][first_input]):
                inputs.append(
                    (
                        f"--labels-file "
                        f"{config['output_dir']}/{wildcards.sample}/norm/"
                        f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                        f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                        f"h5ad2csv/label.csv"
                    )
                )

        elif exists("label_path", config["inputs"][wildcards.sample]):
            inputs.append(
                (
                    f"--labels-file "
                    f"{config['output_dir']}/{wildcards.sample}/norm/"
                    f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                    f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                    f"h5ad2csv/label.csv"
                )
            )

    inputs_str = " ".join(inputs)
    options_str = " ".join(options)
    output_dir = (
        f"{config['output_dir']}/{wildcards.sample}/too-many-cells/{wildcards.maketree}/"
        f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
        f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
    )

    params = (
        f"{inputs_str} {options_str} --output {output_dir} > {output_dir}/cluster.csv"
    )

    return params


def get_maketree_prior_params(wildcards):
    """Parses the options for maketree with a prior.

    Args:
        wildcards: A snakemake object that holds all the wildcards.
    
    Returns:
        A string representing all the options selected for too-many-cells maketree
        using a prior.These options are in the correct format for too-many-cell's 
        command line tool. For example, if this function outputs the string s, then 
        too-many-cells will be executed as such:

        too-many-cells make-tree s

        As an example, s could be:        
        s = "--prior out --labels-file labels.csv --output out2 > clusters.csv"
    """

    # Get the list of user-selected options from the config file
    options = config["too-many-cells"][wildcards.maketree]["options"]

    # Get the first input (input_1)
    first_input = [*config["inputs"].keys()][0]

    if wildcards.norm_method == "tmc":
        if wildcards.sample == "merged":
            if exists("label_path", config["inputs"][first_input]):
                options.insert(
                    0,
                    (
                        f"--labels-file {config['output_dir']}/{wildcards.sample}/"
                        f"too-many-cells/{wildcards.prior}/merged_labels.csv"
                    ),
                )
        elif exists("label_path", config["inputs"][wildcards.sample]):
            options.insert(
                0, f"--labels-file {config['inputs'][wildcards.sample]['label_path']}"
            )

    else:
        if wildcards.sample == "merged":
            if exists("label_path", config["inputs"][first_input]):
                options.insert(
                    0,
                    (
                        f"--labels-file "
                        f"{config['output_dir']}/{wildcards.sample}/norm/"
                        f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                        f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                        f"h5ad2csv/label.csv"
                    ),
                )

        elif exists("label_path", config["inputs"][wildcards.sample]):
            options.insert(
                0,
                (
                    f"--labels-file "
                    f"{config['output_dir']}/{wildcards.sample}/norm/"
                    f"norm_{wildcards.norm_method}_{wildcards.norm_params}/"
                    f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
                    f"h5ad2csv/label.csv"
                ),
            )

    prior = (
         f"{config['output_dir']}/{wildcards.sample}/too-many-cells/{wildcards.prior}/"
         f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
         f"norm_{wildcards.norm_method}_{wildcards.norm_params}"
    )

    options_str = " ".join(options)

    output_dir = (
        f"{config['output_dir']}/{wildcards.sample}/too-many-cells/"
        f"{wildcards.maketree}/"
        f"filter_{wildcards.filter_method}_{wildcards.filter_params}/"
        f"norm_{wildcards.norm_method}_{wildcards.norm_params}"
    )

    params = (
        f"--prior {prior} {options_str} --output {output_dir} "
        f"> {output_dir}/cluster.csv"
    )

    return params


def get_differential_script(wildcards) -> str:
    diff = config['differential'][wildcards.diff]

    if exists("script", diff):
        return diff['script']
    else:
        return "differential.py"


def get_representation(wildcards) -> str:

    if wildcards.c_method == "leiden" or wildcards.c_method == "louvain":
        return "X"
    else:
        return f"X_{wildcards.dr_method}"


