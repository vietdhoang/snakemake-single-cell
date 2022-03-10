# Compute the neighbourhood graph. This is required for some dimensionality
# reduction and clustering methods.
rule neighbourhood_graph:
    input:
        f"{config['output_dir']}/{{sample}}/norm/mtx_{{filter_method}}_{{norm_method}}.h5ad"
    output:
        temp((f"{config['output_dir']}/{{sample}}/dim_reduce/"
              f"mtx_{{filter_method}}_{{norm_method}}.neigh.h5ad"))
    conda:
        "../envs/dim_reduce.yaml"
    script:
        "../scripts/dim_reduce/neighbourhood.py"


rule dim_reduce:
    input:
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
        f"mtx_{{filter_method}}_{{norm_method}}.neigh.h5ad")
    output:
        # dr_method must be a name of a script in the dim_reduce directory
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
        f"mtx_{{filter_method}}_{{norm_method}}_{{dr_method}}.h5ad")
    conda:
        "../envs/dim_reduce.yaml"
    script:
        "../scripts/dim_reduce/{wildcards.dr_method}.py"