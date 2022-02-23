# Compute the neighbourhood graph. This is required for some dimensionality
# reduction and clustering methods.
rule neighbourhood_graph:
    input:
        f"{config['output_dir']}/qc/mtx_{{qc_method}}.h5ad"
    output:
        temp(f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}.neigh.h5ad")
    conda:
        "../envs/dim_reduce.yaml"
    script:
        "../scripts/dim_reduce/neighbourhood.py"


rule dim_reduce:
    input:
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}.neigh.h5ad"
    output:
        # dr_method must be a name of a script in the dim_reduce directory
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    conda:
        "../envs/dim_reduce.yaml"
    script:
        "../scripts/dim_reduce/{wildcards.dr_method}.py"