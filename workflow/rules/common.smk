import os
from typing import List, Union


def make_output_dir(out_dir: Union[str, bytes, os.PathLike]) -> None:
    '''Create an output directory for the pipeline if none exist

    Args:
        out_dir: output path for the directory
    '''
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


def get_final_output() -> List[str]:
    '''Determine the list of desired output files based off what the user
    defined in the config files

    Returns:
        List of strings where each element is an output file
    '''
    
    # Each of these conditions evaluate to True if they exist in the config file
    # and are not empty.
    qc_condition = 'qc_method' in config and config['qc_method']
    dr_condition = 'dim_reduce_method' in config and config['dim_reduce_method']
    cluster_condition = 'cluster_method' in config and config['cluster_method']
    output_condition = 'output' in config and config['output']
    
    final_output = []   # List of desired output files

    # If the user selected qc, dimensionality reduction and clustering
    if qc_condition and dr_condition and cluster_condition:
        
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
    elif qc_condition and dr_condition:

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
    elif qc_condition:

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
    if output_condition:
        final_output.extend(
            expand(
                f"{config['output_dir']}/{{output}}", 
                output=config['other_outputs']
            )
        )
        
    return final_output