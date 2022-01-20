rule mtx_to_h5ad:
    input:
        "{prefix}barcodes.tsv",
        "{prefix}genes.tsv",
        "{prefix}matrix.mtx"
    output:
        config['output_dir'] + "/{prefix}mtx.h5ad"
    conda:
        "../envs/scanpy.yaml"
    params:
        path_in = ".",
        path_out = config['output_dir'] + "/{prefix}mtx.h5ad"
    shell: 
        (   
            f"python {{workflow.basedir}}/scripts/preproc.py mtx_to_h5ad " 
            f"--path={{params.path_in:q}} "
            f"--path_out={{params.path_out:q}} "
            f"--prefix='{{wildcards.prefix}}'"
        )


use rule mtx_to_h5ad as mtx_compressed_to_h5ad with:
    input:
        "{prefix}barcodes.tsv.gz",
        "{prefix}genes.tsv.gz",
        "{prefix}matrix.mtx.gz"


rule qc_seurat:
    input:
        config['output_dir'] + "/{prefix}mtx.h5ad"
    output:
        config['output_dir'] + "/{prefix}mtx_seurat.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (       
            f"python -m fire {{workflow.basedir}}/scripts/qc.py qc_seurat " 
            f"--path='{config['output_dir']}/{{wildcards.prefix}}mtx.h5ad' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}mtx_seurat.h5ad' "
        )