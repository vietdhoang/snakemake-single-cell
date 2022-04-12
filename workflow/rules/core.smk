rule filter:
    input:
        f"{config['output_dir']}/{{sample}}/mtx_orig.h5ad"
    output:
        (f"{config['output_dir']}/{{sample}}/filter/"
         f"filter_{{filter_method}}_{{filter_params}}/mtx.h5ad")
    params:
        script = lambda wildcards: Method('filter', wildcards.filter_method).script,
        params = lambda wildcards: get_params_instance(wildcards.filter_method,
                                                       wildcards.filter_params, 
                                                       get_method_dict('filter'))
    conda:
        "../envs/filter.yaml"
    script:
        "../scripts/filter/{params.script}"


rule norm:
    input:
        (f"{config['output_dir']}/{{sample}}/filter/"
         f"filter_{{filter_method}}_{{filter_params}}/mtx.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/norm/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/mtx.h5ad")
    params:
        script = lambda wildcards: Method('norm', wildcards.norm_method).script,
        params = lambda wildcards: get_params_instance(wildcards.norm_method,
                                                       wildcards.norm_params, 
                                                       get_method_dict('norm'))
    conda:
        "../envs/norm.yaml"
    script:
        "../scripts/norm/{params.script}"


rule dim_reduce:
    input:
        (f"{config['output_dir']}/{{sample}}/norm/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/mtx.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
         f"dim_reduce_{{dr_method}}_{{dr_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/mtx.h5ad")
    params:
        script = lambda wildcards: Method('dim_reduce', wildcards.dr_method).script,
        params = lambda wildcards: get_params_instance(wildcards.dr_method,
                                                       wildcards.dr_params, 
                                                       get_method_dict('dim_reduce'))
    conda:
        "../envs/dim_reduce.yaml"
    script:
        "../scripts/dim_reduce/{params.script}"


rule cluster:
    input:
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
         f"dim_reduce_{{dr_method}}_{{dr_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/mtx.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/cluster/"
         f"cluster_{{c_method}}_{{c_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"dim_reduce_{{dr_method}}_{{dr_params}}/mtx.h5ad")
    params:
        script = lambda wildcards: Method('cluster', wildcards.c_method).script,
        params = lambda wildcards: get_params_instance(wildcards.c_method,
                                                       wildcards.c_params, 
                                                       get_method_dict('cluster')),
        representation = lambda wildcards: f"X_{wildcards.dr_method}"
    conda:
        "../envs/cluster.yaml"
    script:
        "../scripts/cluster/{params.script}"


rule differential:
    input:
        (f"{config['output_dir']}/{{sample}}/cluster/"
         f"cluster_{{c_method}}_{{c_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"dim_reduce_{{dr_method}}_{{dr_params}}/mtx.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/differential/"
         f"differential_{{diff}}/"
         f"filter_{{filter_method}}_{{filter_params}}/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"dim_reduce_{{dr_method}}_{{dr_params}}/"
         f"cluster_{{c_method}}_{{c_params}}/mtx.h5ad")
    params:
        script = "differential/differential.py"
    conda:
        "../envs/differential.yaml"
    script:
        "../scripts/{params.script}"