# The rule for clustering.
# Requires that dimensionality reduction has been run first. 
# Like the dim_reduce rule, it calls a specified function in cluster.py
# to perform the clustering. The function is specified by the c_method wildcard.
rule cluster:
    input:
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
         f"mtx_{{qc_method}}_{{dr_method}}.h5ad")
    output:
        # c_method must be a name of a script in the cluster directory
        (f"{config['output_dir']}/{{sample}}/cluster/"
         f"mtx_{{qc_method}}_{{dr_method}}_{{c_method}}.h5ad")
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


def get_maketree_params(wildcards):
    
    # Get the list of user-selected options from the config file
    options = config['too-many-cells'][wildcards.maketree]['options']
    
    inputs = []

    # If the user decides to not skip the preprocessing of too-many-cells
    if wildcards.qc_method == "tmc_qc":
        # If the sample's name is not "merged"
        if wildcards.sample == "merged":
            for input in config['inputs']:
                inputs.append(f"--matrix-path {config['inputs'][input]['data_path']}")
            
            # Check if the user provided label files by checking
            # if a label file exists for the first input (i.e input_1)
            first_input = [*config['inputs'].keys()][0] 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"--labels-file {config['output_dir']}/merged/"
                     f"too-many-cells/{wildcards.maketree}/merged_labels.csv")
                )
        
        # If the sample's name is not "merged"
        else:
            # Append qc'd matrix csv file to inputs
            inputs.append(
                f"--matrix-path {config['inputs'][wildcards.sample]['data_path']}"
            )

            # Add the label file if provided.
            if exists('label_path', config['inputs'][wildcards.sample]):
                inputs.append(
                    f"--labels-file {config['inputs'][wildcards.sample]['label_path']}"
                )
    
    # If the user decides to use this pipeline's qc methods instead of 
    # too-many-cells' filtering and normalization   
    else: 
        # Append qc'd matrix csv file to inputs
        inputs.append(
            (f"--matrix-path {config['output_dir']}/{wildcards.sample}/qc/"
             f"h5ad2csv/mtx_{wildcards.qc_method}.csv")
        )

        # Append corresponding label file for that csv file (if it exists).
        # If the sample is 'merged' and a label file exists for the first input,
        # then include the --labels-file option for the 'merged' sample
        # If the sample is not 'merged' (i.e anything else), the check if the 
        # sample has a label file. If it does, add it to the list of inputs.
        first_input = [*config['inputs'].keys()][0]
        if (wildcards.sample == "merged"): 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"--labels-file {config['output_dir']}/{wildcards.sample}/"
                     f"qc/h5ad2csv/label_{wildcards.qc_method}.csv")
            )   
        
        elif exists('label_path', config['inputs'][wildcards.sample]):
            inputs.append(
                (f"--labels-file {config['output_dir']}/{wildcards.sample}/qc/"
                 f"h5ad2csv/label_{wildcards.qc_method}.csv")
            )

    inputs_str = " ".join(inputs)
    options_str = " ".join(options)
    output_dir = (f"{config['output_dir']}/{wildcards.sample}/"
                  f"too-many-cells/{wildcards.maketree}/{wildcards.qc_method}")   
    
    params = f"{inputs_str} {options_str} --output {output_dir} > {output_dir}/cluster.csv"

    return params


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

        # Append corresponding label file for that csv file if it exists
        # If the sample is 'merged' and a label file exists for the first input,
        # then it must exist for all other inputs (user must satisfy this)
        # If the sample is not 'merged' (i.e anything else), the check if the 
        # sample has a label file. If it does, add it to the list of inputs.
        first_input = [*config['inputs'].keys()][0]
        if (wildcards.sample == "merged"): 
            if exists('label_path', config['inputs'][first_input]):
                inputs.append(
                    (f"{config['output_dir']}/{wildcards.sample}/qc/"
                     f"h5ad2csv/label_{wildcards.qc_method}.csv")
                )
        
        elif exists('label_path', config['inputs'][wildcards.sample]):
            inputs.append(
                (f"{config['output_dir']}/{wildcards.sample}/qc/"
                f"h5ad2csv/label_{wildcards.qc_method}.csv")
            )

    return inputs


rule tmc_maketree:
    input:
        get_maketree_input
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{maketree}}/{{qc_method}}/{{maketree}}.done")
    
    params:
        parameters = get_maketree_params
    
    shell:
        (      
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}} && "
            f"too-many-cells make-tree {{params.parameters}} && "
            f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}}/{{wildcards.maketree}}.done"
        )


def get_maketree_prior_params(wildcards):
    
    # Get the list of user-selected options from the config file
    options = config['too-many-cells'][wildcards.maketree]['options']
    print(options)
    if wildcards.sample == "merged":
        first_input = [*config['inputs'].keys()][0] 
        if exists('label_path', config['inputs'][first_input]):
            options.insert(0, 
                (f"--labels-file {config['output_dir']}/{wildcards.sample}/"
                 f"qc/h5ad2csv/label_{wildcards.qc_method}.csv")
            )


    elif exists('label_path', config['inputs'][wildcards.sample]):
        options.insert(0,
            (f"--labels-file {config['output_dir']}/{wildcards.sample}/qc/"
             f"h5ad2csv/label_{wildcards.qc_method}.csv")
        )
    
    prior = (f"{config['output_dir']}/{wildcards.sample}/"
             f"too-many-cells/{wildcards.prior}/{wildcards.qc_method}")
    options_str = " ".join(options)
    output_dir = (f"{config['output_dir']}/{wildcards.sample}/"
                  f"too-many-cells/{wildcards.maketree}/{wildcards.qc_method}")

    params = (f"--prior {prior} {options_str} --output {output_dir} "
              f"> {output_dir}/cluster.csv")
    
    return params


rule tmc_maketree_prior:
    input:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{prior}}/{{qc_method}}/{{prior}}.done")        
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{maketree}}/{{qc_method}}/{{prior}}.prior.{{maketree}}.done")
    params:
        parameters = get_maketree_prior_params    
    shell:
        (   
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}} && "   
            f"too-many-cells make-tree {{params.parameters}} && "
            f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}}/{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
        )