rule qc_seurat:
    input:
        config['output_dir'] + "/{prefix}_mtx.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_seurat.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (       
            f"python -m fire {{workflow.basedir}}/scripts/qc.py qc_seurat " 
            f"--path='{config['output_dir']}/{{wildcards.prefix}}_mtx.h5ad' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}_mtx_seurat.h5ad' "
        )