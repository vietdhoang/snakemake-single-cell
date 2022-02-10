# Compute the neighbourhood graph. This is required for some dimensionality
# reduction and clustering methods.
rule neighbourhood_graph:
    input:
        f"{config['output_dir']}/qc/mtx_{{qc_method}}.h5ad"
    output:
        temp(f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}.neigh.h5ad")
    conda:
        "../envs/dim_reduce.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/dim_reduce.py neighbourhood "
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )


# The rule for dimensionality reduction.
# Requires a qc'd file with neighbourhood graph as input, 
# and calls a specified function in dim_reduce.py to perform the dimensionality 
# reduction. The function is specified by the dr_method wildcard.
rule dim_reduce:
    input:
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}.neigh.h5ad"
    output:
        # dr_method must be a name of a function in dim_reduce.py because
        # it will be used to call that function
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    conda:
        "../envs/dim_reduce.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/dim_reduce.py "
            f"{{wildcards.dr_method}} "     # Call the function in dim_reduce.py
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )