# The rule for clustering.
# Requires that dimensionality reduction has been run first. 
# Like the dim_reduce rule, it calls a specified function in cluster.py
# to perform the clustering. The function is specified by the c_method wildcard.
rule cluster:
    input:
        f"{config['output_dir']}/{{sample}}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    output:
        # c_method must be a name of a script in the cluster directory
        f"{config['output_dir']}/{{sample}}/cluster/mtx_{{qc_method}}_{{dr_method}}_{{c_method}}.h5ad"
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/cluster/{wildcards.c_method}.py"


# def get_maketree_params(wildcards):
    
#     options = " ".join(config['too-many-cells'][wildcards.maketree]['options'])
#     if exists('label_file', config):
#         options = f"--labels-file {config['label_file']['file']} {options}"
    
#     output_dir = f"{config['output_dir']}/too-many-cells/{wildcards.maketree}"
#     if hasattr(wildcards, 'prior'):
#         params = (
#             f"--prior {config['output_dir']}/too-many-cells/{wildcards.prior} " 
#             f"{options} "
#             f"--output {output_dir} "
#             f"> {output_dir}/cluster.csv"
#         )
#     else: 
#         params = (
#             f"--matrix-path {config['input']} {options} --output {output_dir} "
#             f"> {output_dir}/cluster.csv"
#         )
        
#     return params


def get_maketree_input(wildcards):

    inputs = []   
    
    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.qc_method == "tmc_qc":
        #  If the sample's name is not "merged"
        if wildcards.sample == "merged":
            # Append the path of all inputs
            # The path to a concatonated label file is required as well.
            for input in config['inputs']:
                inputs.append(config['inputs'][input]['data_path'])
            
            # Check if the user provided label files by checking
            # if a label file exists for the first input (i.e input_1)
            first_input = [*config['inputs'].keys()][0] 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"{config['output_dir']}/merged/too-many-cells/"
                     f"{wildcards.maketree}/merged_labels.csv")
                )
        #  If the sample's name is not "merged" 
        else:
            inputs.append(
                config['inputs'][wildcards.sample]['data_path']
            )

            # Add the label file if provided.
            if exists('label_path', config['inputs'][wildcards.sample]):
                inputs.append(
                    config['inputs'][wildcards.sample]['label_path']
                )

    # If the user decides to use this pipeline's qc methods instead of 
    # too-many-cells' filtering and normalization   
    else: 
        # Append qc'd matrix csv file to inputs
        inputs.append(
            (f"{config['output_dir']}/{wildcards.sample}/qc/"
             f"h5ad2csv/mtx_{wildcards.qc_method}.csv")
        )

        # Append corresponding label file for that csv file
        inputs.append(
            (f"{config['output_dir']}/{wildcards.sample}/qc/"
             f"h5ad2csv/label_{wildcards.qc_method}.csv")
        )

    return inputs


rule tmc_maketree:
    input:
        get_maketree_input
    output:
        f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/{{qc_method}}/{{maketree}}.done"
    
    params:
        # parameters = get_maketree_params
        parameters = ""
    
    shell:
        (      
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/ && "
            f"too-many-cells make-tree {{params.parameters}} && "
            f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.maketree}}.done"
        )


# def get_maketree_prior_input(wildcards):
#     # If the user provided a prior for make-tree, the only input required is 
#     # a flag file indicating that the prior has finished being created.
#     if exists('prior', config['too-many-cells'][wildcards.maketree]):
#         inputs.append(
#             (f"{config['output_dir']}/{wildcards.sample}/too-many-cells/"
#              f"{wildcards.prior}/{wildcards.prior}.done")
#         )

rule tmc_maketree_prior:
    input:
        # get_maketree_input
        f"{config['output_dir']}/{{sample}}/too-many-cells/{{prior}}/{{qc_method}}/{{prior}}.done"        
    output:
        f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/{{qc_method}}/{{prior}}.prior.{{maketree}}.done"
    params:
        # parameters = get_maketree_params
        parameters = ""
    
    shell:
        (   
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/ && "   
            f"too-many-cells make-tree {{params.parameters}} && "
            f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
        )