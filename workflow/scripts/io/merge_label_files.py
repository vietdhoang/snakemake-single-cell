import fire
import pandas as pd
import sys

from os.path import dirname

# Add the scripts directory to Python path and import local files in scripts/
sys.path.insert(0, dirname(dirname(dirname(__file__))))
from scripts.io.merge_samples import modify_obs_name

def merge_label_files(snakemake) -> None:
    
    # Get label list
    label_list = []
    for path_label in list(snakemake.input):
        label_list.append(pd.read_csv(path_label))
    
    # Get list of tags for observations
    obs_tags = []
    for sample in snakemake.config['inputs']:
        obs_tags.append(snakemake.config['inputs'][sample]['obs_tag'])
    
    # Modify observation tags
    for i, df_label in enumerate(label_list):
        df_label['item'] = df_label['item'].map(
            lambda x: modify_obs_name(x, obs_tag=obs_tags[i])
        )
    
    df_merged = pd.concat(label_list)
    df_merged.to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":

    if 'snakemake' in globals():
        merge_label_files(snakemake)
    else:
        fire.Fire()