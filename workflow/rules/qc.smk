# rule quality_control:
#     input:
#         mtx_folder = "data/"
#     output:
#         "figures/filter_genes_dispersion.png",
#         "figures/scatter_mit_gene_expr.png",
#         "figures/scatter_total_gene_count.png",
#         "figures/violin_computed_qc_measures.png"

#     singularity:
#         "docker://vietdhoang/docker-scanpy:1.0"
    
#     shell:
#         "python3 scripts/qc.py"