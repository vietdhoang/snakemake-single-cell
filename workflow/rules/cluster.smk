# The rule for clustering.
# Requires that dimensionality reduction has been run first. 
# Like the dim_reduce rule, it calls a specified function in cluster.py
# to perform the clustering. The function is specified by the c_method wildcard.
rule cluster:
    input:
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    output:
        # c_method must be a name of a function in cluster.py because
        # it will be used to call that function
        f"{config['output_dir']}/cluster/mtx_{{qc_method}}_{{dr_method}}_{{c_method}}.h5ad"
    conda:
        "../envs/cluster.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/cluster.py "
            f"{{wildcards.c_method}} "      # Call the function in cluster.py  
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )


# rule tmc_make_tree:
#     input:
#         config['input']
#     output:
#         (f"{config['output_dir']}/too-many-cells/"
#          f"{config['too-many-cells'][{{mktree_key}}]['output']}/{{mktree_key}}.done")
    
#     params:
#         parameters = get_maketree_params()
    
#     shell:
#         (      
#             f"too-many-cells make-tree {*params.parameters}" 
#             f"> {config[{{wildcards.mktree_key}}]}"
#         )