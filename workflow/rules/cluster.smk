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


# def get_maketree_input(wildcards):
#     [(f"{config['output_dir']}/{wildcards.sample}/"
#                f"{wildcards.pipeline_stage}/{wildcards.basename}/mtx.csv")]


# rule tmc_make_tree:
#     input:
#         get_maketree_input
#     output:
#         f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/{{maketree}}.done"
    
#     params:
#         parameters = get_maketree_params
    
#     shell:
#         (      
#             f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/ && "
#             f"too-many-cells make-tree {{params.parameters}} && "
#             f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.maketree}}.done"
#         )


# rule tmc_make_tree_prior:
#     input:
#         config['input'],
#         f"{config['output_dir']}/{{sample}}/too-many-cells/{{prior}}/{{prior}}.done"
#     output:
#         f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/{{prior}}.prior.{{maketree}}.done"
#     params:
#         parameters = get_maketree_params
    
#     shell:
#         (   
#             f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/ && "   
#             f"too-many-cells make-tree {{params.parameters}} && "
#             f"touch {config['output_dir']}/{{wildcards.sample}}/too-many-cells/{{wildcards.maketree}}/{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
#         )