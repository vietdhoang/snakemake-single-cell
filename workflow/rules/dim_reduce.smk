rule pca:
    input:
        config['output_dir'] + "/{prefix}_mtx_{preproc}.h5ad"
    output:
        config['output_dir'] + "/{prefix}_mtx_{preproc}_pca.h5ad"
    conda:
        "../envs/scanpy.yaml"
    shell:
        (       
            f"python -m fire {{workflow.basedir}}/scripts/dim_reduce.py " 
            f"pca " 
            f"--path='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}.h5ad' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}_pca.h5ad' "
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
            f"python -m fire {{workflow.basedir}}/scripts/dim_reduce.py " 
            f"umap " 
            f"--path='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}.h5ad' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}_umap.h5ad' "
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
            f"python -m fire {{workflow.basedir}}/scripts/dim_reduce.py " 
            f"tsne " 
            f"--path='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}.h5ad' "
            f"--path_out='{config['output_dir']}/{{wildcards.prefix}}_mtx_{{wildcards.preproc}}_tsne.h5ad' "
        )