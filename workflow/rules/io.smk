import os

rule gunzip:    
    input:
        "{prefix}.gz" 
    output:
        f"{config['output_dir']}/{{prefix}}"    
    shell:
        "gunzip -c {input} > {output}"

rule h5ad:
    input:
        config['output_dir'] + "/{prefix}_barcodes.tsv",
        config['output_dir'] + "/{prefix}_genes.tsv",
        config['output_dir'] + "/{prefix}_matrix.mtx"
    output:
        config['output_dir'] + "/{prefix}.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell: 
        (       
            f"python -m fire {{workflow.basedir}}/scripts/io.py mtx_to_h5ad " 
            f"--path='{config['output_dir']}' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}.h5ad' "
            f"--prefix='{{wildcards.prefix}}_'"
        )
        