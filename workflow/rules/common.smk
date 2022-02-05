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


def make_output_dir(dir_out: Union[str, bytes, os.PathLike]) -> None:
    '''Create an output directory for the pipeline if none exist

    Args:
        dir_out: output path for the directory
    '''
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)


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


def get_plots() -> List[str]:
    '''Determine the list of desired plots based off what the user
    defined in the config files

    Returns:
        List of strings where each element is a plot file
    '''

    if not exists('plot_method', config):
        return []

    plot_output = []

    if (exists('scatter', config['plot_method']) 
        and exists('dim_reduce_method', config)
        and exists('dim_reduce', config['plot_method']['scatter'])):
        
        if (exists('label_file', config) 
            and config['plot_method']['scatter']['use_labels']):
            plot_output.extend(
                expand(
                    (f"{config['output_dir']}/figures/labels/"
                     f"{config['output_prefix']}scatter_{{qc}}_{{dr}}_{{l}}.html"),
                    qc=config['qc_method'],
                    dr=config['plot_method']['scatter']['dim_reduce'],
                    l=config['label_file']['labels']
                )                 
            )
            
        else:
            plot_output.extend(
                expand(
                    (f"{config['output_dir']}/figures/no_labels/"
                     f"{config['output_prefix']}scatter_{{qc}}_{{dr}}.html"),
                    qc=config['qc_method'], 
                    dr=config['plot_method']['scatter']['dim_reduce']
                )
            )        
        
        if (exists('cluster', config['plot_method']['scatter']) 
            and exists('cluster_method', config)):            
            plot_output.extend(
                expand(
                    (f"{config['output_dir']}/figures/cluster_assignments/"
                     f"{config['output_prefix']}scatter_{{qc}}_{{dr}}_{{c}}.html"), 
                    qc=config['qc_method'],
                    dr=config['plot_method']['scatter']['dim_reduce'],
                    c=config['plot_method']['scatter']['cluster'] 
                )
            )
    
    return plot_output


def get_final_output() -> List[str]:
    '''Determine the list of desired output files based off what the user
    defined in the config files

    Returns:
        List of strings where each element is an output file
    '''
    
    final_output = []   # List of desired output files

    # If the user selected qc, dimensionality reduction and clustering
    if (exists('qc_method', config) 
        and exists('dim_reduce_method', config) 
        and exists('cluster_method', config)):
        
        # Here we only add the clustering outputs to final_output because
        # Snakemake will automatically produce all other upstream files, such as
        # the dimensionality reduction outputs.
        final_output.extend(
            expand(
                (f"{config['output_dir']}/cluster/"
                 f"{config['output_prefix']}mtx_{{qc}}_{{dimred}}_{{cluster}}.h5ad"),
                qc=config['qc_method'],
                dimred=config['dim_reduce_method'],
                cluster=config['cluster_method']
            )
        )
    
    # If the user selected qc and dimensionality reduction only
    elif (exists('qc_method', config) 
          and exists('dim_reduce_method', config)):

        # Here we only add the dimensionality reduction outputs to final_output
        # for the same reason as above.
        final_output.extend(
            expand(
                (f"{config['output_dir']}/dim_reduce/"
                 f"{config['output_prefix']}mtx_{{qc}}_{{dimred}}.h5ad"),
                qc=config['qc_method'],
                dimred=config['dim_reduce_method']
            )
        )
    
    # If the user selected qc only
    elif (exists('qc_method', config)):

        # Add only qc outputs to final_output
        final_output.extend(
            expand(
                (f"{config['output_dir']}/qc/"
                 f"{config['output_prefix']}mtx_{{qc}}.h5ad"),
                qc=config['qc_method']
            )
        )

    # Takes everything from other_outputs in the config file and appends it
    # to the list of final outputs.
    if (exists('other_outputs', config)):
        final_output.extend(
            expand(
                f"{config['output_dir']}/{{output}}", 
                output=config['other_outputs']
            )
        )
    
    final_output.extend(get_plots())
        
    return final_output


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
