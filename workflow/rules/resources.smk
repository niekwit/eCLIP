rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log",
    cache: False
    threads: 1
    resources:
        runtime=30,
    conda:
        "../envs/mapping.yaml"
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


use rule get_fasta as get_blacklist with:
    output:
        resources.blacklist,
    params:
        url=resources.blacklist_url,
    log:
        "logs/resources/get_blacklist.log",


rule get_yeolab_script:
    output:
        "resources/yeolab/{script}",
    retries: 3
    params:
        url=lambda wildcards: YEOLAB_SCRIPTS[wildcards.script],
    wildcard_constraints:
        script="[A-Za-z0-9_.]+",
    log:
        "logs/resources/get_yeolab_script_{script}.log",
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/mapping.yaml"
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


rule index_fasta:
    input:
        resources.fasta,
    output:
        resources.fai,
    log:
        "logs/resources/index_fasta.log",
    params:
        extra="",  # optional params string
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    wrapper:
        f"{wrapper_version}/bio/samtools/faidx"


rule chrom_sizes:
    input:
        resources.fai,
    output:
        resources.chrom_sizes,
    log:
        "logs/resources/chrom_sizes.log",
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/mapping.yaml"
    shell:
        "cut -f1,2 {input} > {output} 2> {log}"


rule get_repeat_elements:
    # Repeat element reference for the repeat element filtering step of ENCODE eCLIP.
    # ENCODE used RepBase, which is not freely available: this is a free substitute
    # (Dfam consensus sequences + rDNA repeating unit)
    output:
        resources.repeat_fasta,
    retries: 3
    params:
        clade=resources.dfam_clade,
        rdna=resources.rdna_accession,
    log:
        "logs/resources/get_repeat_elements.log",
    cache: False
    threads: 1
    resources:
        runtime=30,
    conda:
        "../envs/mapping.yaml"
    script:
        "../scripts/get_repeat_elements.py"


rule star_index_genome:
    input:
        fasta=resources.fasta,
        fai=resources.fai,
        gtf=resources.gtf,
    output:
        directory(resources.star_index),
    params:
        gtf=lambda wildcards, input: input.gtf,
    log:
        "logs/resources/star_index_genome.log",
    cache: False
    threads: config["resources"]["star_index"]["cpu"]
    resources:
        runtime=config["resources"]["star_index"]["time"],
        mem_mb=40000,
    conda:
        "../envs/mapping.yaml"
    script:
        "../scripts/star_index.sh"


rule star_index_repeats:
    input:
        fasta=resources.repeat_fasta,
        fai=f"{resources.repeat_fasta}.fai",
    output:
        directory(resources.star_repeat_index),
    params:
        gtf="",  # no annotation
    log:
        "logs/resources/star_index_repeats.log",
    cache: False
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["star_index"]["time"],
        mem_mb=8000,
    conda:
        "../envs/mapping.yaml"
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
    params:
        fasta=adapters.barcodes_fasta(),
    log:
        "logs/resources/barcodes_fasta.log",
    localrule: True
    conda:
        "../envs/mapping.yaml"
    shell:
        "printf '%s' '{params.fasta}' > {output} 2> {log}"
