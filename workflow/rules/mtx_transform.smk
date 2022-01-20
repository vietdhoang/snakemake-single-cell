wildcard_constraints:
    # A prefix can be anything that doesn't contain a '\', '/', or a whitespace
    prefix = "[^\\/\s]*"

rule pca:
    input:
        config['output_dir'] + "/{prefix}_mtx_{preproc}.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_pca.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/mtx_transform.py pca "  
            f"--path={{input:q}} "
            f"--path_out={{output:q}} "
        )


rule umap:
    input:
        config['output_dir'] + "/{prefix}_mtx_{preproc}.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_umap.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/mtx_transform.py umap "  
            f"--path={{input:q}} "
            f"--path_out={{output:q}} "
        )


rule tsne:
    input:
        config['output_dir'] + "/{prefix}_mtx_{preproc}.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_tsne.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/mtx_transform.py tsne "  
            f"--path={{input:q}} "
            f"--path_out={{output:q}} "
        )

rule leiden:
    input:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_umap.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_leiden.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/mtx_transform.py leiden "  
            f"--path={{input:q}} "
            f"--path_out={{output:q}} "
        )