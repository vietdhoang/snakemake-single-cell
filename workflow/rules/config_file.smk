import itertools
import os
from typing import Final, List, Union


def get_core_output() -> List[str]:
    '''Determine the list of desired output based off what the user defined in
    the config files. The core output includes output from filtering,
    normalization, dimensionality reduction and clustering.

    Returns:
        List of strings where each element is an output file from furtherest
        point downstream in the pipeline. For example, if the user only 
        wants dimensionality reduction, then the output list will contain 
        outputs of dimensionality reduction. If the user wants clustering,
        the the output list will contain outputs for clustering instead since
        clustering is further downstream. There is no need to add the output
        for dimensionality reduction since that will automatically be generated
        anyway.
    '''

    core_output = []

    if exists('cluster_method', config):
        # Here we only add the clustering outputs to final_output because
        # Snakemake will automatically produce all other upstream files, such as
        # the dimensionality reduction outputs.
        core_output.extend(
            expand(
                ("cluster/mtx_{filter}_{norm}_{dimred}_{cluster}.h5ad"),
                filter=config['filter_method'],
                norm=config['norm_method'],
                dimred=config['dim_reduce_method'],
                cluster=config['cluster_method']
            )
        )
    
    elif exists('dim_reduce_method', config):
        # Here we only add the dimensionality reduction outputs to final_output
        # for the same reason as above.
        core_output.extend(
            expand(
                ("dim_reduce/mtx_{filter}_{norm}_{dimred}.h5ad"),
                filter=config['filter_method'],
                norm=config['norm_method'],
                dimred=config['dim_reduce_method']
            )
        )
    
    elif exists('norm_method', config):
        core_output.extend(
            expand(
                ("norm/mtx_{filter}_{norm}.h5ad"),
                filter=config['filter_method'],
                norm=config['norm_method']
            )
        )
    
    elif exists('filter_method', config):
        core_output.extend(
            expand(
                ("filter/mtx_{filter}_{norm}.h5ad"),
                filter=config['filter_method']
            )
        )

    return core_output
    

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
                    "figures/labels/scatter_{filter}_{norm}_{dr}_{l}.html",
                    filter=config['filter_method'],
                    norm=config['norm_method'],
                    dr=config['plot_method']['scatter']['dim_reduce'],
                    l=config['label_file']['labels']
                )                 
            )
            
        else:
            plot_output.extend(
                expand(
                    "figures/no_labels/scatter_{filter}_{norm}_{dr}.html",
                    filter=config['filter_method'],
                    norm=config['norm_method'],
                    dr=config['plot_method']['scatter']['dim_reduce']
                )
            )        
        
        if (exists('cluster', config['plot_method']['scatter']) 
            and exists('cluster_method', config)):            
            plot_output.extend(
                expand(
                    "figures/cluster_assignments/scatter_{filter}_{norm}_{dr}_{c}.html", 
                    filter=config['filter_method'],
                    norm=config['norm_method'],
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

    samples = [*config['inputs'].keys()]
    if len(samples) > 1 and config['run_on_each_sample']:
        samples.append('merged')
    elif len(samples) > 1 and not config['run_on_each_sample']:
        samples = ['merged']

    for key in dict_tmc:
        if 'make-tree' in key:            
            # If the make-tree uses a prior
            if exists('prior', dict_tmc[key]):                
                if dict_tmc[dict_tmc[key]['prior']]['tmc_qc']:
                    tmc_output.extend(
                        expand(
                            (f"{config['output_dir']}/{{sample}}/too-many-cells/"
                             f"{key}/tmc_filter/tmc_norm/"
                             f"{dict_tmc[key]['prior']}.prior.{key}.done"),
                            sample=samples
                        )
                    )
                else:
                    tmc_output.extend(
                        expand(
                            (f"{config['output_dir']}/{{sample}}/too-many-cells/"
                             f"{key}/{{filter_method}}/{{norm_method}}"
                             f"{dict_tmc[key]['prior']}.prior.{key}.done"),
                            sample=samples,
                            filter_method=config['filter_method'],
                            norm_method=config['norm_method'] 
                        )
                    )
                
            # If the make-tree doesn't contain a prior
            else:
                if dict_tmc[key]['tmc_qc']:
                    tmc_output.extend(
                        expand(
                            (f"{config['output_dir']}/{{sample}}/too-many-cells/"
                             f"{key}/tmc_filter/tmc_norm/{key}.done"),
                            sample=samples
                        )
                    )
                else:
                    tmc_output.extend(
                        expand(
                            (f"{config['output_dir']}/{{sample}}/too-many-cells/"
                             f"{key}/{{filter_method}}/{{norm_method}}/{key}.done"),
                            sample=samples,
                            filter_method = config['filter_method'],
                            norm_method = config['norm_method']
                        )
                    )
    
    return tmc_output


def get_final_output() -> List[str]:
    '''Determine the list of desired output files based off what the user
    defined in the config files

    Returns:
        List of strings where each element is an output file
    '''
    
    final_output = []   # List of desired output files

    # Get core output (filtering, normalization, dim reduction, clustering)
    final_output.extend(get_core_output())    

    # Get output plots
    final_output.extend(get_plots())

    # Takes everything from other_outputs in the config file and appends it
    # to the list of final outputs.
    if (exists('other_outputs', config)):
        final_output.extend(
            expand("{output}", output=config['other_outputs'])
        )

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

    # Get too-many-cells output
    final_output.extend(get_too_many_cells_output())
        
    return final_output