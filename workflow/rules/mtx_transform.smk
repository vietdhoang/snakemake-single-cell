# The rule for dimensionality reduction.
# Requires a qc'd file as input, and calls a specified function in dim_reduce.py
# to perform the dimensionality reduction. The function is specified by the
# dr_method wildcard.
rule dim_reduce:
    input:
        config['output_dir'] + "/qc/{prefix}mtx_{qc_method}.h5ad"
    output:
        # dr_method must be a name of a function in dim_reduce.py because
        # it will be used to call that function
        config['output_dir'] + "/dim_reduce/{prefix}mtx_{qc_method}_{dr_method}.h5ad"
    conda:
        "../envs/dim_reduce.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/dim_reduce.py "
            f"{{wildcards.dr_method}} "     # Call the function in dim_reduce.py
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )


# The rule for clustering.
# Requires that dimensionality reduction has been run first. 
# Like the dim_reduce rule, it calls a specified function in cluster.py
# to perform the clustering. The function is specified by the c_method wildcard.
rule cluster:
    input:
        config['output_dir'] + "/dim_reduce/{prefix}mtx_{qc_method}_{dr_method}.h5ad"
    output:
        # c_method must be a name of a function in cluster.py because
        # it will be used to call that function
        config['output_dir'] + "/cluster/{prefix}mtx_{qc_method}_{dr_method}_{c_method}.h5ad"
    conda:
        "../envs/cluster.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/cluster.py "
            f"{{wildcards.c_method}} "      # Call the function in cluster.py  
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )