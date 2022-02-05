ruleorder: mtx_to_h5ad_genes > mtx_to_h5ad 


# Converts gene-barcode matricies to h5ad format
rule mtx_to_h5ad:
    input:
        config['input']    
    output:
        config['output_dir'] + "/{prefix}mtx.h5ad"
    conda:
        "../envs/preproc.yaml"
    params:
        path_in = ".",
        path_out = config['output_dir'] + "/{prefix}mtx.h5ad",
        path_label = config['label_file']['file'] if exists('label_file', config) else None,
        labels = list_to_str(config['label_file']['labels']) if exists('label_file', config) else None
    shell: 
        (   
            f"python {{workflow.basedir}}/scripts/preproc.py mtx_to_h5ad " 
            f"--path={{params.path_in:q}} "
            f"--path_out={{params.path_out:q}} "
            f"--prefix={{wildcards.prefix:q}} "
            f"--path_label={{params.path_label:q}} "
            f"--labels={{params.labels}}"
        )


# Adapts mtx_to_h5ad to work with uncompressed gene-barcode matrices
# This is part of a work-around for a bug on scanpy where it cannot process
# compressed gene-barcode matrices. This rule comes right after gunzip.
use rule mtx_to_h5ad as mtx_to_h5ad_genes with:
    input:
        "{prefix}barcodes.tsv",
        "{prefix}genes.tsv",
        "{prefix}matrix.mtx"


# Uncompresses gene-barcode matricies
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


# The rule for quality control.
# Requires an h5ad file as input, and calls a specified function in qc.py
# to perform the quality control. The function is specified by the
# qc_method wildcard.
rule qc:
    input:
        config['output_dir'] + "/{prefix}mtx.h5ad"
    output:
        # qc_method must be a name of a function in qc.py because
        # it will be used to call that function
        config['output_dir'] + "/qc/{prefix}mtx_{qc_method}.h5ad"
    conda:
        "../envs/preproc.yaml"
    shell:
        (      
            f"python {{workflow.basedir}}/scripts/preproc.py "
            f"{{wildcards.qc_method}} "     # Call the function in qc.py 
            f"--path_in={{input:q}} "
            f"--path_out={{output:q}} "
        )