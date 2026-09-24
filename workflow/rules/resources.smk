rule get_fasta:
    output:
        resources.fasta,
    log:
        "logs/resources/get_fasta.log",
    cache: False
    retries: 3
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        runtime=30,
    params:
        url=resources.fasta_url,
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
    output:
        resources.gtf,
    log:
        "logs/resources/get_gtf.log",
    params:
        url=resources.gtf_url,


use rule get_fasta as get_blacklist with:
    output:
        resources.blacklist,
    log:
        "logs/resources/get_blacklist.log",
    params:
        url=resources.blacklist_url,


rule get_yeolab_script:
    output:
        "resources/yeolab/{script}",
    log:
        "logs/resources/get_yeolab_script_{script}.log",
    wildcard_constraints:
        script="[A-Za-z0-9_.]+",
    retries: 3
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        url=lambda wildcards: YEOLAB_SCRIPTS[wildcards.script],
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


rule index_fasta:
    input:
        resources.fasta,
    output:
        resources.fai,
    log:
        "logs/resources/index_fasta.log",
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    params:
        extra="",  # optional params string
    wrapper:
        f"{wrapper_version}/bio/samtools/faidx"


rule chrom_sizes:
    input:
        resources.fai,
    output:
        resources.chrom_sizes,
    log:
        "logs/resources/chrom_sizes.log",
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        runtime=10,
    shell:
        "cut -f1,2 {input} > {output} 2> {log}"


rule get_repeat_elements:
    # Repeat element reference for the repeat element filtering step of ENCODE eCLIP.
    # ENCODE used RepBase, which is not freely available: this is a free substitute
    # (Dfam consensus sequences + rDNA repeating unit)
    output:
        resources.repeat_fasta,
    log:
        "logs/resources/get_repeat_elements.log",
    cache: False
    retries: 3
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        runtime=30,
    params:
        clade=resources.dfam_clade,
        rdna=resources.rdna_accession,
    script:
        "../scripts/get_repeat_elements.py"


rule star_index_genome:
    input:
        fasta=resources.fasta,
        fai=resources.fai,
        gtf=resources.gtf,
    output:
        directory(resources.star_index),
    log:
        "logs/resources/star_index_genome.log",
    cache: False
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["star_index"]["cpu"]
    resources:
        runtime=config["resources"]["star_index"]["time"],
        mem_mb=40000,
    params:
        gtf=lambda wildcards, input: input.gtf,
    script:
        "../scripts/star_index.sh"


rule star_index_repeats:
    input:
        fasta=resources.repeat_fasta,
        fai=f"{resources.repeat_fasta}.fai",
    output:
        directory(resources.star_repeat_index),
    log:
        "logs/resources/star_index_repeats.log",
    cache: False
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["star_index"]["time"],
        mem_mb=8000,
    params:
        gtf="",  # no annotation
    script:
        "../scripts/star_index.sh"


use rule index_fasta as index_repeats with:
    input:
        resources.repeat_fasta,
    output:
        f"{resources.repeat_fasta}.fai",
    log:
        "logs/resources/index_repeats.log",


rule barcodes_fasta:
    # Barcodes for demultiplexing of paired-end reads by eclipdemux
    output:
        "resources/yeolab/barcodes.fasta",
    log:
        "logs/resources/barcodes_fasta.log",
    localrule: True
    conda:
        "../envs/mapping.yaml"
    params:
        fasta=adapters.barcodes_fasta(),
    shell:
        "printf '%s' '{params.fasta}' > {output} 2> {log}"
