rule tmc_maketree:
    input:
        get_maketree_input
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/{{maketree}}.done")    
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
            f"filter_{{wildcards.filter_method}}_{{wildcards.filter_params}}/"
            f"norm_{{wildcards.norm_method}}_{{wildcards.norm_params}}/ && "

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"filter_{{wildcards.filter_method}}_{{wildcards.filter_params}}/"
            f"norm_{{wildcards.norm_method}}_{{wildcards.norm_params}}/"
            f"{{wildcards.maketree}}.done"
        )


rule tmc_maketree_prior:
    input:    
        (f"{config['output_dir']}/{{sample}}/too-many-cells/{{prior}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/{{prior}}.done")        
    output:
        (f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"{{prior}}.prior.{{maketree}}.done")
    params:
        parameters = get_maketree_prior_params
    shell:
        (   
            f"mkdir -p {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"filter_{{wildcards.filter_method}}_{{wildcards.filter_params}}/"
            f"norm_{{wildcards.norm_method}}_{{wildcards.norm_params}}/ && "            

            f"too-many-cells make-tree {{params.parameters}} && "

            f"touch {config['output_dir']}/{{wildcards.sample}}/"
            f"too-many-cells/{{wildcards.maketree}}/"
            f"filter_{{wildcards.filter_method}}_{{wildcards.filter_params}}/"
            f"norm_{{wildcards.norm_method}}_{{wildcards.norm_params}}/"
            f"{{wildcards.prior}}.prior.{{wildcards.maketree}}.done"
        )