import os

rule gunzip:    
    input:
        in_bar = "{prefix}_barcodes.tsv.gz",
        in_gene = "{prefix}_genes.tsv.gz",
        in_mat = "{prefix}_matrix.mtx.gz"    
    output:
        out_bar = f"{config['output_dir']}/{{prefix}}_barcodes.tsv",
        out_gene = f"{config['output_dir']}/{{prefix}}_genes.tsv",
        out_mat = f"{config['output_dir']}/{{prefix}}_matrix.mtx"    
    shell:
        "gunzip -c {input.in_bar} > {output.out_bar}"
        "&& gunzip -c {input.in_gene} > {output.out_gene}"
        "&& gunzip -c {input.in_mat} > {output.out_mat}"


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
        