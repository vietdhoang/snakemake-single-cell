rule plot_scatter_labels:
    input:
        f"{config['output_dir']}/dim_reduce/mtx_{{qc_method}}_{{dr_method}}.h5ad"
    output:
        expand(
            (f"{config['output_dir']}/figures/labels/"
             "scatter_{{qc_method}}_{{dr_method}}_{label}.html"),
            label=config['label_file']['labels']
        )
    conda:
        "../envs/plot.yaml"
    params:
        path_dir_out = f"{config['output_dir']}/figures/labels/",
        prefix = f"scatter_{{qc_method}}_{{dr_method}}",
        labels = list_to_str(config['label_file']['labels']),
        use_rep = f"{{dr_method}}"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/plot.py "
            f"scatter "
            f"--path_in={{input:q}} "
            f"--path_dir_out={{params.path_dir_out:q}} "
            f"--prefix={{params.prefix:q}} "
            f"--labels={{params.labels}} "
            f"--use_rep={{params.use_rep:q}}"
        )


use rule plot_scatter_labels as plot_scatter_no_labels with:
    output:
        (f"{config['output_dir']}/figures/no_labels/"
         f"scatter_{{qc_method}}_{{dr_method}}.html")
    params:
        path_dir_out = f"{config['output_dir']}/figures/no_labels/",
        prefix = f"scatter_{{qc_method}}_{{dr_method}}",
        labels = list_to_str([]),
        use_rep = f"{{dr_method}}"


use rule plot_scatter_labels as plot_scatter_clusters with:
    input:
        (f"{config['output_dir']}/cluster/"
         f"mtx_{{qc_method}}_{{dr_method}}_{{c_method}}.h5ad")
    output:
        (f"{config['output_dir']}/figures/cluster_assignments/"
         f"scatter_{{qc_method}}_{{dr_method}}_{{c_method}}.html")
    params:    
        path_dir_out = f"{config['output_dir']}/figures/cluster_assignments/",
        prefix = f"scatter_{{qc_method}}_{{dr_method}}",
        labels = list_to_str(["{c_method}"]),
        use_rep = f"{{dr_method}}"