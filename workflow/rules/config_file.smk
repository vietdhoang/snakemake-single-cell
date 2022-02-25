import itertools
import os
from typing import Final, List, Union

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
                    "figures/labels/scatter_{qc}_{dr}_{l}.html",
                    qc=config['qc_method'],
                    dr=config['plot_method']['scatter']['dim_reduce'],
                    l=config['label_file']['labels']
                )                 
            )
            
        else:
            plot_output.extend(
                expand(
                    "figures/no_labels/scatter_{qc}_{dr}.html",
                    qc=config['qc_method'], 
                    dr=config['plot_method']['scatter']['dim_reduce']
                )
            )        
        
        if (exists('cluster', config['plot_method']['scatter']) 
            and exists('cluster_method', config)):            
            plot_output.extend(
                expand(
                    "figures/cluster_assignments/scatter_{qc}_{dr}_{c}.html", 
                    qc=config['qc_method'],
                    dr=config['plot_method']['scatter']['dim_reduce'],
                    c=config['plot_method']['scatter']['cluster'] 
                )
            )
    
    return plot_output


def get_too_many_cells_output() -> List[str]:
    '''Determine the list of desired plots based off what the user
    defined in the config files

    Returns:
        List of strings where each element is a plot file
    '''

    if not exists('too-many-cells', config):
        return []

    tmc_output = []
    dict_tmc = config['too-many-cells']

    samples = [*config[inputs].keys()]
    if len(sample) > 1 and config['run_on_each_sample']:
        samples.append('merged')
    elif len(sample) > 1 and not config['run_on_each_sample']:
        samples = ['merged']

    for key in dict_tmc:
        if 'make-tree' in key:            
            # If only prior exists but not comb_options
            if exists('prior', dict_tmc[key]):                
                tmc_output.append(
                    expand(
                        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
                         f"{key}/{dict_tmc[key]['prior']}.prior.{key}.done"),
                        sample=samples
                    )
                )
            # If neither comb_options nor prior exist
            else:
                tmc_output.append(
                    (f"{config['output_dir']}/too-many-cells/"
                     f"{key}/{key}.done")
                )
    
    return tmc_output


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
                ("cluster/mtx_{qc}_{dimred}_{cluster}.h5ad"),
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
                "dim_reduce/mtx_{qc}_{dimred}.h5ad",
                qc=config['qc_method'],
                dimred=config['dim_reduce_method']
            )
        )
    
    # If the user selected qc only
    elif (exists('qc_method', config)):

        # Add only qc outputs to final_output
        final_output.extend(
            expand(
                "qc/mtx_{qc}.h5ad",
                qc=config['qc_method']
            )
        )

    # Takes everything from other_outputs in the config file and appends it
    # to the list of final outputs.
    if (exists('other_outputs', config)):
        final_output.extend(
            expand("{output}", output=config['other_outputs'])
        )
    

    final_output.extend(get_plots())


    # Prepend the output directory path and sample name to each of the 
    # elements in the list of final outputs
    samples = [*config['inputs'].keys()]
        
    if len(samples) > 1:    
        if config['run_on_each_sample']:
            samples.append('merged')
            final_output = [
                f"{config['output_dir']}/{sample}/{file}" 
                for sample in samples 
                for file in final_output
            ]
        
        else:
            final_output= [
                f"{config['output_dir']}/merged/{file}" 
                for file in final_output
            ]
    
    else:
        final_output = [
            f"{config['output_dir']}/{samples[0]}/{file}" 
            for file in final_output
        ]

    final_output.extend(get_too_many_cells_output())
        
    return final_output