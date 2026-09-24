import os
import re
import glob
import datetime
import itertools
import pandas as pd
from scripts.resources import GenomeResources
import scripts.adapters as adapters
from snakemake.utils import min_version, validate
from snakemake.logging import logger

# Yeo lab scripts of the ENCODE pipeline that are not available via conda,
# these are downloaded at pinned commits (https://github.com/YeoLab/eclip, https://github.com/YeoLab/merge_peaks)
_ECLIP = "https://raw.githubusercontent.com/YeoLab/eclip/c0fffc4979a92371dc0667a03e3d957bf7f77600/bin"
_MERGE_PEAKS = "https://raw.githubusercontent.com/YeoLab/merge_peaks/aedc0a14d4ba109ee65678a3201a52c5bb6ad473/bin"

YEOLAB_SCRIPTS = {
    "barcodecollapsepe.py": f"{_ECLIP}/barcodecollapsepe.py",
    "overlap_peakfi_with_bam.pl": f"{_ECLIP}/overlap_peakfi_with_bam.pl",
    "compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl": f"{_ECLIP}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl",
    "compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl": f"{_MERGE_PEAKS}/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl",
    "make_informationcontent_from_peaks.pl": f"{_MERGE_PEAKS}/perl/make_informationcontent_from_peaks.pl",
    "parse_idr_peaks.pl": f"{_MERGE_PEAKS}/perl/parse_idr_peaks.pl",
    "get_reproducing_peaks.pl": f"{_MERGE_PEAKS}/perl/get_reproducing_peaks.pl",
    "full_to_bed.py": f"{_MERGE_PEAKS}/full_to_bed.py",
}


def yeolab_script(name):
    """
    Returns path of downloaded Yeo lab script
    """
    return f"resources/yeolab/{name}"


def paired_end():
    """
    Checks if paired-end or single-end reads are used
    """
    sample = csv["sample"].tolist()[0]

    if os.path.isfile(f"reads/{sample}_R1_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    elif os.path.isfile(f"reads/{sample}.fastq.gz"):
        logger.info("Single-end reads detected...")
        return False
    else:
        raise ValueError(
            f"No read files found for sample {sample}, expected reads/{sample}.fastq.gz (single-end) "
            f"or reads/{sample}_R1_001.fastq.gz and reads/{sample}_R2_001.fastq.gz (paired-end)"
        )


def samples(PAIRED_END):
    """
    Checks samples.csv and returns IP samples, input samples and all samples
    """
    ALL = csv["sample"].tolist()
    IP = csv.loc[csv["control"] != "", "sample"].tolist()
    INPUT = list(dict.fromkeys(csv.loc[csv["control"] != "", "control"].tolist()))

    # Sample names must be unique
    duplicated = sorted(set(x for x in ALL if ALL.count(x) > 1))
    if len(duplicated) != 0:
        raise ValueError("Following samples are not unique:\n" + "\n".join(duplicated))

    # IP samples: end with _[0-9]+, and can not contain _vs_ (used to name IDR comparisons)
    wrong = [x for x in IP if not re.match(r"^.+_[0-9]+$", x)]
    if len(wrong) != 0:
        raise ValueError(
            "Following IP samples do not end with _[replicate number]:\n"
            + "\n".join(wrong)
        )
    wrong = [x for x in ALL if "_vs_" in x]
    if len(wrong) != 0:
        raise ValueError("Sample names can not contain _vs_:\n" + "\n".join(wrong))

    # All input samples must have their own row without a control
    missing = [x for x in INPUT if x not in ALL]
    if len(missing) != 0:
        raise ValueError(
            "Following input samples are not listed in the sample column:\n"
            + "\n".join(missing)
        )
    wrong = [x for x in INPUT if x in IP]
    if len(wrong) != 0:
        raise ValueError(
            "Following samples are used as input, but have a control themselves:\n"
            + "\n".join(wrong)
        )
    unused = [x for x in ALL if x not in IP and x not in INPUT]
    if len(unused) != 0:
        raise ValueError(
            "Following samples have no control and are not used as input for another sample:\n"
            + "\n".join(unused)
        )

    # Check if sample names match file names
    not_found = []
    for sample in ALL:
        if PAIRED_END:
            files = [
                f"reads/{sample}_R1_001.fastq.gz",
                f"reads/{sample}_R2_001.fastq.gz",
            ]
        else:
            files = [f"reads/{sample}.fastq.gz"]
        not_found.extend([f for f in files if not os.path.isfile(f)])
    if len(not_found) != 0:
        raise ValueError("Following files not found:\n" + "\n".join(not_found))

    # Paired-end reads need inline barcode IDs
    if PAIRED_END:
        for col in ["barcode_a", "barcode_b"]:
            if col not in csv.columns:
                raise ValueError(
                    f"Paired-end reads require column {col} in samples.csv"
                )
        wrong = csv.loc[
            (csv["barcode_a"] == "") | (csv["barcode_b"] == ""), "sample"
        ].tolist()
        if len(wrong) != 0:
            raise ValueError(
                "Following samples have no barcode IDs (use NIL for barcode-less libraries):\n"
                + "\n".join(wrong)
            )

    return IP, INPUT, ALL


def sample_barcodes(sample):
    """
    Returns unique inline barcode IDs of a paired-end sample
    """
    row = csv.loc[csv["sample"] == sample].iloc[0]

    return list(dict.fromkeys([row["barcode_a"], row["barcode_b"]]))


def units():
    """
    Returns a dictionary with the library units (wildcard unit) and their sample.
    A unit is a sample (single-end) or a sample/barcode combination (paired-end),
    as paired-end reads are demultiplexed by their inline barcode.
    """
    _units = {}
    for sample in SAMPLES:
        if PAIRED_END:
            for barcode in sample_barcodes(sample):
                _units[f"{sample}.{barcode}"] = sample
        else:
            _units[sample] = sample

    return _units


def condition(sample):
    """
    Returns condition of IP sample (sample name without replicate number)
    """
    return re.sub(r"_[0-9]+$", "", sample)


def idr_pairs():
    """
    Returns list of IP replicate pairs (of the same condition) for IDR analysis
    """
    pairs = []
    conditions = list(dict.fromkeys(condition(x) for x in IP_SAMPLES))
    for cond in conditions:
        replicates = [x for x in IP_SAMPLES if condition(x) == cond]
        if len(replicates) < 2:
            logger.info(
                f"WARNING: condition {cond} has only one replicate: no IDR analysis for this condition..."
            )
        pairs.extend(itertools.combinations(replicates, 2))

    return pairs


def control(sample):
    """
    Returns control sample (size-matched input) of IP sample
    """
    return csv.loc[csv["sample"] == sample, "control"].iloc[0]


def se_adapter(sample):
    """
    Returns adapter set of single-end sample
    """
    if "adapter" in csv.columns:
        adapter = csv.loc[csv["sample"] == sample, "adapter"].iloc[0]
        if adapter:
            return adapter

    return config["cutadapt"]["se_adapter"]


def cutadapt_args(unit, round_):
    """
    Returns adapter arguments of cutadapt as string for a unit.

    Round 1 removes both 5' and 3' adapters, round 2 only the 3' adapters of read 2
    (single-end: same 3' adapters as in round 1) to control for double ligation events.
    """
    sample = UNITS[unit]

    if PAIRED_END:
        barcode = unit.split(".")[-1]
        barcodes = sample_barcodes(sample)
        args = []
        if round_ == 1:
            args.append(f"-a {adapters.PE_R1_ADAPTER}")
            args.extend(f"-g {x}" for x in adapters.pe_read1_5p_adapters(barcodes))
        args.extend(f"-A {x}" for x in adapters.pe_read2_adapters(barcodes))
    else:
        args = [f"-a {x}" for x in adapters.se_adapters(se_adapter(sample))]

    return " ".join(args)


def star_input(wildcards):
    """
    Returns read files for STAR mapping of a unit
    """
    return {
        "reads": expand(
            "results/repeats/unmapped/{unit}.r{end}.fq.gz",
            unit=wildcards.unit,
            end=ENDS,
        )
    }


def peaks_input(wildcards):
    """
    Returns IP and input BAM files and read counts for input normalisation
    """
    input_sample = control(wildcards.sample)

    return {
        "ip_bam": f"results/mapped/{wildcards.sample}.bam",
        "ip_bai": f"results/mapped/{wildcards.sample}.bam.bai",
        "ip_num": f"results/mapped/{wildcards.sample}.readnum.txt",
        "input_bam": f"results/mapped/{input_sample}.bam",
        "input_bai": f"results/mapped/{input_sample}.bam.bai",
        "input_num": f"results/mapped/{input_sample}.readnum.txt",
    }


def multiqc_input():
    """
    Returns input files for MultiQC
    """
    files = []

    # FastQC of raw reads
    if PAIRED_END:
        files.extend(
            expand(
                "results/qc/fastqc/raw/{sample}_R{end}_fastqc.zip",
                sample=SAMPLES,
                end=ENDS,
            )
        )
    else:
        files.extend(
            expand("results/qc/fastqc/raw/{sample}_fastqc.zip", sample=SAMPLES)
        )

    # FastQC after each cutadapt round
    files.extend(
        expand(
            "results/qc/fastqc/trim{round}/{unit}.r{end}_fastqc.zip",
            round=[1, 2],
            unit=UNITS,
            end=ENDS,
        )
    )

    # cutadapt reports and STAR logs
    files.extend(
        expand("results/qc/cutadapt/round{round}/{unit}.txt", unit=UNITS, round=[1, 2])
    )
    files.extend(expand("results/star/repeats/{unit}.Log.final.out", unit=UNITS))
    files.extend(expand("results/star/genome/{unit}.Log.final.out", unit=UNITS))

    return files


def targets():
    """
    Returns file targets for rule all
    """
    TARGETS = ["results/qc/multiqc.html"]

    # Bigwig files of all libraries (IP and input)
    TARGETS.extend(
        expand(
            "results/bigwig/{sample}.norm.{strand}.bw",
            sample=SAMPLES,
            strand=["pos", "neg"],
        )
    )

    # Input normalised peaks of each IP sample
    TARGETS.extend(
        expand(
            "results/peaks/{sample}/{sample}.peaks.{ext}",
            sample=IP_SAMPLES,
            ext=["bed", "narrowPeak", "bb"],
        )
    )
    TARGETS.extend(
        expand("results/peaks/{sample}/{sample}.total_entropy.txt", sample=IP_SAMPLES)
    )

    # Reproducible peaks of replicate pairs
    TARGETS.extend(
        [
            f"results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.reproducible_peaks.{ext}"
            for s1, s2 in IDR_PAIRS
            for ext in ["bed", "custombed", "narrowPeak", "bb"]
        ]
    )

    return TARGETS
