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
            (f"{config['output_dir']}/input_{i}/qc/"
            f"mtx_{wildcards.qc_method}.h5ad")
        )
        i += 1
            
    return inputs


def parse_tmc_options_other(dict_make_tree: dict) -> List[str]:
    return " ".join(dict_make_tree['other_options'])


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


# def permute_comb_vals(dict_make_tree: dict) -> List[dict]:
#     keys = [*dict_make_tree['comb_options'].keys()]
#     values = [*dict_make_tree['comb_options'].values()]

#     return [dict(zip(keys, v)) for v in itertools.product(*values)]