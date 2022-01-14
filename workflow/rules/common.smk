import os
import sys
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
    
    # Get the output directory and prefix from config file
    output_dir = config['output_dir']
    output_prefix = config['output_prefix']
    
    final_output = []
    
    # Loop through the outputs.
    for output in config['outputs']:

        # If the output contains all, then go to the 'all' entry in 
        # all_outputs.yaml and add all possible output files.
        if output == 'all':
            final_output.extend(
                expand(f"{output_dir}/{output_prefix}{{dataset}}", 
                       dataset=config['all'])
            )

            return final_output        

    
    final_output.extend(
        expand(f"{output_dir}/{output_prefix}{{dataset}}", 
               dataset=config['outputs'])
    )
    
    return final_output