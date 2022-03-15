ruleorder: merge_samples > normalize

# Converts gene-barcode matricies to h5ad format
rule mtx_to_h5ad:
    input:
        lambda wildcards: config['inputs'][wildcards.sample]['data_path']    
    output:
        f"{config['output_dir']}/{{sample}}/mtx.h5ad"
    conda:
        "../envs/preproc.yaml"
    script: 
        "../scripts/io/mtx_to_h5ad.py"


rule h5ad_to_csv:
    input:
        f"{config['output_dir']}/{{sample}}/{{pipeline_stage}}/mtx_{{basename}}.h5ad"
    output:
        (f"{config['output_dir']}/{{sample}}/"
         f"{{pipeline_stage}}/h5ad2csv/mtx_{{basename}}.csv"),
        (f"{config['output_dir']}/{{sample}}/"
         f"{{pipeline_stage}}/h5ad2csv/label_{{basename}}.csv")
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/io/h5ad_to_csv.py"


rule filter:
    input:
        f"{config['output_dir']}/{{sample}}/mtx.h5ad"
    output:
        f"{config['output_dir']}/{{sample}}/filter/mtx_{{filter_method}}.h5ad"
    params:
        filter_method = lambda wildcards: wildcards.filter_method
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/filter/filter.py"


rule normalize:
    input:
        f"{config['output_dir']}/{{sample}}/filter/mtx_{{filter_method}}.h5ad"
    output:
        # norm_method must be a name of a script in the norm directory
        f"{config['output_dir']}/{{sample}}/norm/mtx_{{filter_method}}_{{norm_method}}.h5ad"
    params:
        norm_method = wildcards.norm_method,
        log = True
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/normalize/normalize.py"


rule merge_samples:
    input:
        get_merge_samples_input
    output:
        f"{config['output_dir']}/merged/norm/mtx_{{filter_method}}_{{norm_method}}.h5ad"
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/io/merge_samples.py"


rule merge_label_files:
    input:
        lambda wildcards: [config['inputs'][sample]['label_path'] for sample in [*config['inputs'].keys()]]
    output:
        f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/merged_labels.csv"
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/io/merge_label_files.py"



