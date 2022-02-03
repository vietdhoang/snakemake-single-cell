import os
from typing import List, Union


def make_output_dir(dir_out: Union[str, bytes, os.PathLike]) -> None:
    '''Create an output directory for the pipeline if none exist

    Args:
        dir_out: output path for the directory
    '''
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)


def exists(key: str, dictionary: dict) -> bool:
    '''Helper function to determine if a key exists in a 
    dictionary and isn't empty
    '''
    return key in dictionary and dictionary[key]


def get_plots() -> List[str]:
    '''Determine the list of desired plots based off what the user
    defined in the config files

    Returns:
        List of strings where each element is a plot file
    '''

    plot_output = []

    if (exists('scatter', config['plot_method']) 
        and exists('dim_reduce', config)
        and exists('dim_reduce', config['plot_method']['scatter'])):
        
        if exists('label_file', config):
            for dr in config['plot_method']['scatter']['dim_reduce']:
                for label in config['label_file']['label']:
                    plot_output.append(
                        (f"{config['output_dir']}/figures/labels"
                         f"{config['output_prefix']}scatter_{dr}_{label}.html")
                    )
            
        else:
            plot_output.extend(
                expand(
                    (f"{config['output_dir']}/figures/dim_reduce"
                     f"{config['output_prefix']}scatter_{{dr}}.html"), 
                    dr=configconfig['plot_method']['scatter']['dim_reduce']
                )
            )        
        
        if (exists('cluster', config['plot_method']['scatter']) 
            and exists('cluster', config)):            
            plot_output.extend(
                expand(
                    (f"{config['output_dir']}/figures/cluster_assignments"
                     f"{config['output_prefix']}scatter_{{dr}}_{{c}}.html"), 
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
          and exists('dim_reduce_method', config):

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
     elif (exists('qc_method', config):

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
    if (exists('other_outputs', config):
        final_output.extend(
            expand(
                f"{config['output_dir']}/{{output}}", 
                output=config['other_outputs']
            )
        )
    
    final_output.extend(get_plots())
        
    return final_output