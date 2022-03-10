import os
from typing import List, Union


def get_wc_constraint(key: str) -> str:
    '''Determine wildcard constraints for wildcards that are dependent on
    entries in the config file. For example, the qc wildcard should be
    constrained to the entries under qc_method in the config file.

    Args:
        key: A key for the config dictionary
    
    Returns:
        A regex string that defines the wild constraint associated with
        config[key]
    '''    
    if exists(key, config):
        constraint = "|".join(config[key])
        return f"({constraint})"
    
    # If key does not exist in the config file, return the regex below instead.
    # This is the default regex: anything that doesn't contain a 
    # '\', '.',  '/', or a whitespace
    else:
        return "[^\\/\.\s]*"


def exists(key: str, dictionary: dict) -> bool:
    '''Helper function to determine if a key exists in a 
    dictionary and isn't empty.

    Args:
        key: A key that is potentially in the dictionary
        dictionary: the dictionary that will be searched
    
    Returns:
        True if the key exists in the dictionary and the entry is not empty.
    '''
    return key in dictionary and dictionary[key]


def list_to_str(l: List[str]) -> str:
    '''Helper function that converts a list to a string literal. This is used
    to help Fire to parse lists as input in the command line. See rule scatter
    in rules/plot.smk for an example.

    Args:
        l: list of strings that will be converted into a string
    
    Returns:
        The list as a string
    '''
    l_str = ",".join(l)
    return f"\"[{l_str}]\""


def get_merge_samples_input(wildcards):
    inputs = []
    i = 1
    while exists(f"input_{i}", config['inputs']):
        inputs.append(
            (f"{config['output_dir']}/input_{i}/norm/"
            f"mtx_{wildcards.filter_method}_{wildcards.norm_method}.h5ad")
        )
        i += 1
            
    return inputs


def get_h5ad_to_csv_output(wildcards):

    output = [(f"{config['output_dir']}/{wildcards.sample}/"
               f"{wildcards.pipeline_stage}/{wildcards.basename}_mtx.csv")]

    output.extend(
        expand(
            (f"{config['output_dir']}/{wildcards.sample}/"
             f"{wildcards.pipeline_stage}/{wildcards.basename}_{{label}}.csv"),
            label = config['labels']
        )
    )

    return output


def get_maketree_input(wildcards):

    inputs = []   
    
    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.norm_method == "tmc_norm":
        #  If the sample's name is not "merged"
        if wildcards.sample == "merged":
            # Append the path of all inputs
            # The path to a concatonated label file is required as well.
            for input in config['inputs']:
                inputs.append(config['inputs'][input]['data_path'])
            
            # Check if the user provided label files by checking
            # if a label file exists for the first input (i.e input_1)
            first_input = [*config['inputs'].keys()][0] 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"{config['output_dir']}/merged/too-many-cells/"
                     f"{wildcards.maketree}/merged_labels.csv")
                )
        #  If the sample's name is not "merged" 
        else:
            inputs.append(
                config['inputs'][wildcards.sample]['data_path']
            )

            # Add the label file if provided.
            if exists('label_path', config['inputs'][wildcards.sample]):
                inputs.append(
                    config['inputs'][wildcards.sample]['label_path']
                )

    # If the user decides to use this pipeline's qc methods instead of 
    # too-many-cells' filtering and normalization   
    else: 
        # Append qc'd matrix csv file to inputs
        inputs.append(
            (f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
             f"mtx_{wildcards.filter_method}_{wildcards.norm_method}.csv")
        )

        # Append corresponding label file for that csv file if it exists
        # If the sample is 'merged' and a label file exists for the first input,
        # then it must exist for all other inputs (user must satisfy this)
        # If the sample is not 'merged' (i.e anything else), the check if the 
        # sample has a label file. If it does, add it to the list of inputs.
        first_input = [*config['inputs'].keys()][0]
        if (wildcards.sample == "merged"): 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                     f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv")
                )
        
        elif exists('label_path', config['inputs'][wildcards.sample]):
            inputs.append(
                (f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv")
            )

    return inputs


def get_maketree_params(wildcards):
    
    # Get the list of user-selected options from the config file
    options = config['too-many-cells'][wildcards.maketree]['options']
    
    inputs = []

    # Get the first input (input_1)
    first_input = [*config['inputs'].keys()][0]

    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.norm_method == "tmc_norm":
        # If the sample's name is not "merged"
        if wildcards.sample == "merged":
            for input in config['inputs']:
                inputs.append(f"--matrix-path {config['inputs'][input]['data_path']}")
            
            # Check if the user provided label files by checking
            # if a label file exists for the first input             
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"--labels-file {config['output_dir']}/merged/"
                     f"too-many-cells/{wildcards.maketree}/merged_labels.csv")
                )
        
        # If the sample's name is not "merged"
        else:
            # Append qc'd matrix csv file to inputs
            inputs.append(
                f"--matrix-path {config['inputs'][wildcards.sample]['data_path']}"
            )

            # Add the label file if provided.
            if exists('label_path', config['inputs'][wildcards.sample]):
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
                f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                f"mtx_{wildcards.filter_method}_{wildcards.norm_method}.csv"
            )
        )

        # Append corresponding label file for that csv file (if it exists).
        # If the sample is 'merged' and a label file exists for the first input,
        # then include the --labels-file option for the 'merged' sample
        # If the sample is not 'merged' (i.e anything else), the check if the 
        # sample has a label file. If it does, add it to the list of inputs.
        if (wildcards.sample == "merged"): 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (
                        f"--labels-file "
                        f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                        f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv"
                    )
            )   
        
        elif exists('label_path', config['inputs'][wildcards.sample]):
            inputs.append(
                (
                    f"--labels-file "
                    f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                    f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv"
                )
            )

    inputs_str = " ".join(inputs)
    options_str = " ".join(options)
    output_dir = (
        f"{config['output_dir']}/{wildcards.sample}/too-many-cells/"
        f"{wildcards.maketree}/{wildcards.filter_method}/{wildcards.norm_method}"
    )   
    
    params = f"{inputs_str} {options_str} --output {output_dir} > {output_dir}/cluster.csv"

    return params


def get_maketree_prior_params(wildcards):

    # Get the list of user-selected options from the config file
    options = config['too-many-cells'][wildcards.maketree]['options']
    
    # Get the first input (input_1)
    first_input = [*config['inputs'].keys()][0]
    
    if wildcards.norm_method == "tmc_norm":    
        if wildcards.sample == "merged": 
            if exists('label_path', config['inputs'][first_input]):
                options.insert(0, 
                    (
                        f"--labels-file {config['output_dir']}/{wildcards.sample}/"
                        f"too-many-cells/{wildcards.prior}/merged_labels.csv"
                    )
                )
        elif exists('label_path', config['inputs'][wildcards.sample]):
            options.insert(0,
                f"--labels-file {config['inputs'][wildcards.sample]['label_path']}"
            )

    else:
        if (wildcards.sample == "merged"): 
            if exists('label_path', config['inputs'][first_input]):
                options.insert(0,
                    (
                        f"--labels-file "
                        f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                        f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv"
                    )
                )
            
        elif exists('label_path', config['inputs'][wildcards.sample]):
            options.insert(0,
                (
                    f"--labels-file "
                    f"{config['output_dir']}/{wildcards.sample}/norm/h5ad2csv/"
                    f"label_{wildcards.filter_method}_{wildcards.norm_method}.csv")
            )
    
    prior = (
        f"{config['output_dir']}/{wildcards.sample}/too-many-cells/"
        f"{wildcards.prior}/{wildcards.filter_method}/{wildcards.norm_method}"
    )
    
    options_str = " ".join(options)
    
    output_dir = (
        f"{config['output_dir']}/{wildcards.sample}/too-many-cells/"
        f"{wildcards.maketree}/{wildcards.filter_method}/{wildcards.norm_method}"
    )

    params = (f"--prior {prior} {options_str} --output {output_dir} "
              f"> {output_dir}/cluster.csv")
    
    return params