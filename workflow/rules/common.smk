import os
import sys
from typing import List


def make_output_dir(out_dir: str) -> None:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


def get_final_output() -> List[str]:
    
    output_dir = config['output_dir']
    output_prefix = config['output_prefix']
    
    final_output = []
    
    for output in config['outputs']:
        if output == 'all':
            final_output.extend(
                expand(f"{output_dir}/{output_prefix}{{dataset}}", 
                       dataset=config['all'])
            )

            return final_output        

    
    final_output.extend(
        expand(f"{output_dir}/{output_prefix}{{dataset}}", 
               dataset=config['all'])
    )
    
    return final_output