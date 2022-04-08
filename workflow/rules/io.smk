ruleorder: merge_samples > norm

# Converts gene-barcode matricies to h5ad format
rule mtx_to_h5ad:
    input:
        lambda wildcards: config['inputs'][wildcards.sample]['data_path']    
    output:
        f"{config['output_dir']}/{{sample}}/mtx_orig.h5ad"
    conda:
        "../envs/io.yaml"
    script: 
        "../scripts/io/mtx_to_h5ad.py"


rule h5ad_to_csv:
    input:
        f"{config['output_dir']}/{{sample}}/{{pipeline_stage}}/{{basename}}/mtx.h5ad"
    output:
        (f"{config['output_dir']}/{{sample}}/"
         f"{{pipeline_stage}}/{{basename}}/h5ad2csv/mtx.csv"),
        (f"{config['output_dir']}/{{sample}}/"
         f"{{pipeline_stage}}/{{basename}}/h5ad2csv/label.csv")
    wildcard_constraints:
        pipeline_stage = "filter|norm|cluster|dim_reduce|too-many-cells"
    conda:
        "../envs/io.yaml"
    script:
        "../scripts/io/h5ad_to_csv.py"


rule merge_samples:
    input:
        get_merge_samples_input
    output:
        (f"{config['output_dir']}/merged/norm/"
         f"norm_{{norm_method}}_{{norm_params}}/"
         f"filter_{{filter_method}}_{{filter_params}}/mtx.h5ad")
    conda:
        "../envs/io.yaml"
    script:
        "../scripts/io/merge_samples.py"


rule merge_label_files:
    input:
        lambda wildcards: [config['inputs'][sample]['label_path'] 
                           for sample in [*config['inputs'].keys()]]
    output:
        f"{config['output_dir']}/{{sample}}/too-many-cells/{{maketree}}/merged_labels.csv"
    conda:
        "../envs/io.yaml"
    script:
        "../scripts/io/merge_label_files.py"