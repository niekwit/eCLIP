# since I cannot get HOMER annotation to work for non TSS peaks, I will use this approach: https://www.biostars.org/p/329227/

rule convert_gtf2bed: #only contains transcript lines
    input:
        gtf=resources.gtf
    output:
        bed=resources.bed
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotation/gtf2bed.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        r"""
        gtf2bed {params.extra} < {input.gtf} | awk -F '\t' '{{if ($8~"transcript") print $0}}' > {output.bed} 2> {log}
        """
    

rule annotate_crosslink_sites:
    input:
        bed="results/pureclip/crosslink_sites/{ip_sample}_vs_{input_sample}.bed",
        genome_bed=resources.bed,
    output:
        bed="results/annotation/{ip_sample}_vs_{input_sample}.annotated_crosslinks.bed",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotation/{ip_sample}_vs_{input_sample}.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        "bedtools intersect "
        "-a {input.bed} "
        "-b {input.genome_bed} "
        "-wb "
        "> {output.bed} 2> {log}"


rule extract_geneIDs:
    input:
        bed="results/annotation/{ip_sample}_vs_{input_sample}.annotated_crosslinks.bed",
    output:
        txt="results/annotation/{ip_sample}_vs_{input_sample}.ensembl_gene_ids.txt",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotation/{ip_sample}_vs_{input_sample}.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        r"""
        awk -F '\t' '{{print $11}}' {input.bed} | sort | uniq > {output.txt} 2> {log}
        """
    
    