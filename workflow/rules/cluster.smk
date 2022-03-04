localrules: tmc_maketree_prior

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


rule tmc_maketree:
    input:
        get_maketree_input
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{maketree}}/{{qc_method}}/{{maketree}}.done")
    
    params:
        parameters = get_maketree_params
    resources:
        partition = "veryhimem",
        mem_mb = 100000,
        time = "0-05:00:00"    
    shell:
        (      
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}} && "

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}}/"
            f"{{wildcards.maketree}}.done"
        )


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
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}} && "

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/{{wildcards.qc_method}}/"
            f"{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
        )