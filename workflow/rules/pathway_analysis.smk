rule convert_ensembl_geneID2symbol:
    input:
        "results/annotation/{ip_sample}_vs_{input_sample}.ensembl_gene_ids.txt",
    output:
        "results/annotation/{ip_sample}_vs_{input_sample}.gene_symbols.txt",
    params:
        dataset=resources.ensembl,
        # use 'asia.ensembl.org', 'useast.ensembl.org' if default fails
        host="www.ensembl.org", 
        mart="ENSEMBL_MART_ENSEMBL",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/gseapy/biomart/{ip_sample}_vs_{input_sample}.convert_ensembl2symbol.log",
    conda:
        "../envs/annotation.yaml",
    shell:
        "gseapy biomart "
        "-f ensembl_gene_id {input} "
        "-a external_gene_name "
        "-o {output} "
        "-d {params.dataset} "
        "--host {params.host} "
        "-m {params.mart} "
        "--verbose "
        "> {log} 2>&1"


rule pathway_analysis:
    input:
        "results/annotation/{ip_sample}_vs_{input_sample}.gene_symbols.txt",
    output:
        directory("results/pathway_analysis/{ip_sample}_vs_{input_sample}/"),
    params:
        libs="GO_Biological_Process_2023,GO_Molecular_Function_2023,KEGG_2021_Human,Reactome_2022",
        genome=resources.enrichr,
        fdr=config["fdr_cutoff"],
        terms=10,
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/gseapy/enrichr/{ip_sample}_vs_{input_sample}.pathway_analysis.log",
    conda:
        "../envs/annotation.yaml",
    shell:
        "gseapy enrichr "
        "-i {input} "
        "-g {params.libs} "
        "--organism {params.genome} "
        "--description {wildcards.ip_sample}_vs_{wildcards.input_sample} "
        "--cut-off {params.fdr} "
        "--top-term {params.terms} "
        "--outdir {output} "
        "{params.extra} "
        "--verbose "
        "> {log} 2>&1"


    