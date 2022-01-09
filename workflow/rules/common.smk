import os
import sys
from typing import List, NoReturn


def make_output_dir(out_dir: str) -> NoReturn:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


def get_final_output() -> List[str]:
    final_output = []

    if 'feature-barcode-matrix' in config:
        final_output.extend(
            expand(f"{config['output_dir']}/{{dataset}}", 
                   dataset=config['feature-barcode-matrix'])
        )
    
    if config['h5ad_file']:
        final_output.append(f"{config['output_dir']}/{config['h5ad_file']}")

    
    return final_output