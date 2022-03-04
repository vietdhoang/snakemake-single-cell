localrules: plot_scatter_labels, plot_scatter_no_labels, plot_scatter_clusters

rule plot_scatter_labels:
    input:
        f"{config['output_dir']}/{{sample}}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    output:
        expand(
            (f"{config['output_dir']}/"
             "{{sample}}/figures/labels/scatter_{{qc_method}}_{{dr_method}}_{label}.html"),
            label=config['labels']
        )
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot/scatter.py"


use rule plot_scatter_labels as plot_scatter_no_labels with:
    output:
        (f"{config['output_dir']}/{{sample}}/figures/no_labels/"
         f"scatter_{{qc_method}}_{{dr_method}}.html")


use rule plot_scatter_labels as plot_scatter_clusters with:
    input:
        (f"{config['output_dir']}/{{sample}}/cluster/"
         f"mtx_{{qc_method}}_{{dr_method}}_{{c_method}}.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/figures/cluster_assignments/"
         f"scatter_{{qc_method}}_{{dr_method}}_{{c_method}}.html")