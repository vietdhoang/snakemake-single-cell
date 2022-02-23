# Converts gene-barcode matricies to h5ad format
rule mtx_to_h5ad:
    input:
        config['input']    
    output:
        f"{config['output_dir']}/mtx.h5ad"
    conda:
        "../envs/preproc.yaml"
    script: 
        "../scripts/io.py"

rule qc:
    input:
        f"{config['output_dir']}/mtx.h5ad"
    output:
        # qc_method must be a name of a script in the qc directory
        f"{config['output_dir']}/qc/mtx_{{qc_method}}.h5ad"
    conda:
        "../envs/preproc.yaml"
    script:
        "../scripts/qc/{wildcards.qc_method}.py"
