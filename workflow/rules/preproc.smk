localrules: mtx_to_h5ad, h5ad_to_csv, merge_samples, merge_label_files
ruleorder: merge_samples > qc

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


rule qc:
    input:
        f"{config['output_dir']}/{{sample}}/mtx.h5ad"
    output:
        # qc_method must be a name of a script in the qc directory
        f"{config['output_dir']}/{{sample}}/qc/mtx_{{qc_method}}.h5ad"
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/qc/{wildcards.qc_method}.py"


rule merge_samples:
    input:
        get_merge_samples_input
    output:
        f"{config['output_dir']}/merged/qc/mtx_{{qc_method}}.h5ad"
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



