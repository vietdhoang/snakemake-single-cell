rule plot_scatter_labels:
    input:
        (f"{config['output_dir']}/{{sample}}/dim_reduce/"
         f"mtx_{{filter_method}}_{{norm_method}}_{{dr_method}}.h5ad")
    output:
        expand(
            (f"{config['output_dir']}/"
             "{{sample}}/figures/labels/"
             "scatter_{{filter_method}}_{{norm_method}}_{{dr_method}}_{label}.html"),
            label=config['labels']
        )
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot/scatter.py"


use rule plot_scatter_labels as plot_scatter_dim_reduce with:
    output:
        (f"{config['output_dir']}/{{sample}}/figures/dim_reduce/"
         f"scatter_{{filter_method}}_{{norm_method}}_{{dr_method}}.html")


use rule plot_scatter_labels as plot_scatter_cluster with:
    input:
        (f"{config['output_dir']}/{{sample}}/cluster/"
         f"mtx_{{filter_method}}_{{norm_method}}_{{dr_method}}_{{c_method}}.h5ad")
    output:
        (f"{config['output_dir']}/{{sample}}/figures/cluster/"
         f"scatter_{{filter_method}}_{{norm_method}}_{{dr_method}}_{{c_method}}.html")


# rule plot_after_filter:
#     input:
#         f"{config['output_dir']}/{{sample}}/filter/mtx_{{filter_method}}.h5ad"
#     output:
#         expand(
#             (f"{config['output_dir']}/"
#              "{{sample}}/after_filter/{{filter_method}}/{plot}.html"),
#             plot=["qc_metrics_cells", 
#                   "qc_metrics_genes", 
#                   "scatter_MT_vs_total_count",
#                   "scatter_n_genes_by_count_vs_total_count"]
#         )
#     params:
#         dir_out = lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/after_filter/{wildcards.filter_method}/"     
#     conda:
#         "../envs/plot.yaml"
#     script:
#         "../scripts/plot/qc_metrics.py"


rule plot_after_filter:
    input:
        f"{config['output_dir']}/{{sample}}/filter/mtx_{{filter_method}}.h5ad"
    output:
        expand(
            (f"{config['output_dir']}/"
             "{{sample}}/after_filter/{{filter_method}}/{plot}.html"),
            plot=["density_cells_n_genes_by_counts", 
                  "density_cells_total_counts", 
                  "density_cells_pct_counts_mt",
                  "density_genes_n_cells_by_counts"
                  "density_genes_total_counts"]
        )
    params:
        dir_out = lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/after_filter/{wildcards.filter_method}/"     
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot/qc_metrics.py"


use rule plot_after_filter as plot_before_filter with:
    input:
        f"{config['output_dir']}/{{sample}}/mtx.h5ad"
    output:
        expand(
            (f"{config['output_dir']}/"
             "{{sample}}/before_filter/{{filter_method}}/{plot}.html"),
             plot=["density_cells_n_genes_by_counts", 
                  "density_cells_total_counts", 
                  "density_cells_pct_counts_mt",
                  "density_genes_n_cells_by_counts"
                  "density_genes_total_counts"]
        )
    params:
        dir_out = lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/before_filter/{wildcards.filter_method}/"