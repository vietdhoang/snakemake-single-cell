rule tmc_maketree:
    input:
        get_maketree_input
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{maketree}}/{{filter_method}}/{{norm_method}}/{{maketree}}.done")
    
    params:
        parameters = get_maketree_params
    resources:
        partition = "veryhimem",
        mem_mb = 100000,
        time = "0-05:00:00"    
    shell:
        (      
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"{{wildcards.filter_method}}/{{wildcards.norm_method}} && "

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/{{wildcards.filter_method}}/"
            f"{{wildcards.norm_method}}/{{wildcards.maketree}}.done"
        )


rule tmc_maketree_prior:
    input:    
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{prior}}/{{filter_method}}/{{norm_method}}/{{prior}}.done")        
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/"
         f"{{maketree}}/{{filter_method}}/{{norm_method}}/"
         f"{{prior}}.prior.{{maketree}}.done")
    params:
        parameters = get_maketree_prior_params
    shell:
        (   
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"{{wildcards.filter_method}}/{{wildcards.norm_method}} && "

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"{{wildcards.filter_method}}/{{wildcards.norm_method}}/"
            f"{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
        )