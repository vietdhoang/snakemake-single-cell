# Converts gene-barcode matricies to h5ad format
rule mtx_to_h5ad:
    input:
        config['input']    
    output:
        f"{config['output_dir']}/mtx.h5ad"
    conda:
        "../envs/preproc.yaml"
    params:
        path_label = config['label_file']['file'] if exists('label_file', config) else None,
        labels = list_to_str(config['label_file']['labels']) if exists('label_file', config) else None
    shell: 
        (   
            f"python {{workflow.basedir}}/scripts/preproc.py mtx_to_h5ad " 
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
            f"--path_label={{params.path_label:q}} "
            f"--labels={{params.labels}}"
        )


# The rule for quality control.
# Requires an h5ad file as input, and calls a specified function in qc.py
# to perform the quality control. The function is specified by the
# qc_method wildcard.
rule qc:
    input:
        f"{config['output_dir']}/mtx.h5ad"
    output:
        # qc_method must be a name of a function in qc.py because
        # it will be used to call that function
        f"{config['output_dir']}/qc/mtx_{{qc_method}}.h5ad"
    conda:
        "../envs/preproc.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/preproc.py "
            f"{{wildcards.qc_method}} "     # Call the function in qc.py 
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )