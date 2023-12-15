rule get_fasta:
    output:
        resources.fasta
    retries: 3 # retry up to 3 times in case of network failure
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    conda:
        "../envs/trim.yaml",
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz > {log} 2>&1"


rule get_gtf: 
    output:
        resources.gtf
    retries: 3 # retry up to 3 times in case of network failure
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log"
    conda:
        "../envs/trim.yaml",
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz > {log} 2>&1"


rule install_homer_genome:
    output:
        touch("resources/homer_genome_installed"),
    params:
        genome=config["genome"],
    log:
        "logs/resources/homer_install_genome.log"
    conda:
        "../envs/pureclip.yaml"
    script:
        "../scripts/install_homer_genome.sh"
