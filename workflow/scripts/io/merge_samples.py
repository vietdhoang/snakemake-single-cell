import fire
import os
import pandas as pd
import scanpy as sc

from typing import Union


def modify_obs_name(obs_name: str, obs_tag: str = "") -> str:

    if not obs_tag:
        return obs_name
    elif '-' in obs_name:
        return f"{obs_name[:obs_name.find('-')]}-{obs_tag}"
    else:
        return f"{obs_name}-{obs_tag}"


def get_sample_name_from_path(path: Union[str, bytes, os.PathLike], 
                              config: dict) -> str:
    sample_list = [*config['inputs']]
    for sample_name in sample_list:
        if sample_name in path:
            return sample_name
    
    return ""


def merge_samples(snakemake) -> None:
    '''Merge two samples into one
    '''
    adatas = []
    obs_tags = []
    
    for path_sample in list(snakemake.input):        
        adata = sc.read_h5ad(path_sample)
        sample_name = get_sample_name_from_path(path_sample, snakemake.config)
        
        obs_names = pd.Series(adata.obs_names, name='item')
        obs_tag = snakemake.config['inputs'][sample_name]['obs_tag']
        adata.obs_names = obs_names.map(lambda x: modify_obs_name(x, obs_tag=obs_tag))
        
        adatas.append(adata)
        obs_tags.append(obs_tag)

    merged_adata = adatas[0].concatenate(
        *adatas[1:], 
        join='outer',
        batch_categories=obs_tags, 
        fill_value=0,
        index_unique=None
    )
    
    merged_adata.write(snakemake.output[0])   


if __name__ == "__main__":

    if 'snakemake' in globals():
        merge_samples(snakemake)    
    else:
        fire.Fire()