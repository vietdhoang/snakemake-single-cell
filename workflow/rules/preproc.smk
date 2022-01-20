wildcard_constraints:
    # A prefix can be anything that doesn't contain a '\', '/', or a whitespace
    prefix = "[^\\/\s]*"


rule gunzip:    
    input:
        in_bar = "{prefix}barcodes.tsv.gz",
        in_gene ="{prefix}genes.tsv.gz",
        in_mtx = "{prefix}matrix.mtx.gz" 
    output:
        out_bar = temp("{prefix}barcodes.tsv"),
        out_gene = temp("{prefix}genes.tsv"),
        out_mtx = temp("{prefix}matrix.mtx")    
    shell:
        (
            "gunzip -c {input.in_bar} > {output.out_bar} "
            "&& gunzip -c {input.in_gene} > {output.out_gene} "
            "&& gunzip -c {input.in_mtx} > {output.out_mtx} "
        )


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


use rule mtx_to_h5ad as mtx_to_h5ad_features_compressed with:
    input:
        "{prefix}barcodes.tsv.gz",
        "{prefix}features.tsv.gz",
        "{prefix}matrix.mtx.gz"


use rule mtx_to_h5ad as mtx_to_h5ad_features_uncompressed with:
    input:
        "{prefix}barcodes.tsv",
        "{prefix}features.tsv",
        "{prefix}matrix.mtx"


rule qc_seurat:
    input:
        config['output_dir'] + "/{prefix}mtx.h5ad"
    output:
        config['output_dir'] + "/{prefix}mtx_qcseurat.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (       
            f"python {{workflow.basedir}}/scripts/preproc.py qc_seurat " 
            f"--path={{input:q}} "
            f"--path_out={{output:q}}"
        )