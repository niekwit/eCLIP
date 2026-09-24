"""
Adapter and barcode sequences used by the ENCODE eCLIP processing pipeline
(eCLIP-seq Processing Pipeline v2.2, https://www.encodeproject.org/pipelines/ENCPL357ADL/)

Cutadapt is run with 15 nt "chunks" of each adapter, so that adapters with truncated
5' ends (which cutadapt would not find by default) are also trimmed. This follows
https://github.com/YeoLab/eclip (generate_adaptertrim_fasta.ipynb)
"""

# Illumina adapter sequences (3' adapters of read 1 and read 2)
R1_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R2_ADAPTER = "AGATCGGAAGAGCGTCGTGT"

# 3' adapter of read 1 in paired-end eCLIP (as in the ENCODE SOP; the 5 Ns match the random-mer)
PE_R1_ADAPTER = "NNNNN" + R1_ADAPTER

# 5' adapter of read 1 in paired-end eCLIP (barcode is ligated to the 3' end of this sequence)
PE_5P_ADAPTER = "CTTCCGATCT"

# Inline barcodes of paired-end eCLIP (yeolabbarcodes_20170101.fasta). NIL = barcode-less (size-matched input)
PE_BARCODES = {
    "A01": "AAGCAAT",
    "B06": "GGCTTGT",
    "C01": "ACAAGTT",
    "D8f": "TGGTCCT",
    "A03": "ATGACCNNNNT",
    "G07": "TCCTGTNNNNT",
    "A04": "CAGCTTNNNNT",
    "F05": "GGATACNNNNT",
    "NIL": "",
}

# 3' adapters of single-end eCLIP libraries: barcode-like prefix + Illumina read 1 adapter
# (InvRNA1-8_adapters.fasta and InvRil19_adapters.yaml)
SE_ADAPTER_PREFIXES = {
    "InvRil19": "",
    "InvRNA1": "AGCGCTAG",
    "InvRNA2": "GATATCGA",
    "InvRNA3": "CGCAGACG",
    "InvRNA4": "TATGAGTA",
    "InvRNA5": "AGGTGCGT",
    "InvRNA6": "GAACATAC",
    "InvRNA7": "ACATAGCG",
    "InvRNA8": "GTGCGATA",
}

CHUNK_SIZE = 15


def revcomp(seq):
    """Returns the reverse complement of a sequence (IUPAC N allowed)"""
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def sliding_windows(seq, size=CHUNK_SIZE):
    """Returns all subsequences of a given length (5' to 3')"""
    return [seq[i : i + size] for i in range(len(seq) - size + 1)]


def se_adapters(name):
    """
    Returns 3' adapter chunks for single-end eCLIP reads (cutadapt -a)

    The adapters that have a barcode-like prefix are padded with Ns at their 5' end
    for the first two chunks
    """
    if name not in SE_ADAPTER_PREFIXES:
        raise ValueError(
            f"Unknown single-end adapter {name}, choose from: {', '.join(SE_ADAPTER_PREFIXES)}"
        )

    prefix = SE_ADAPTER_PREFIXES[name]
    full = prefix + R1_ADAPTER

    chunks = []
    if prefix:
        chunks.extend(
            [
                "NN" + full[: CHUNK_SIZE - 2],
                "N" + full[: CHUNK_SIZE - 1],
            ]
        )
    chunks.extend(sliding_windows(full))

    return chunks


def pe_read2_adapters(barcode_ids):
    """
    Returns 3' adapter chunks of read 2 for paired-end eCLIP (cutadapt -A)

    These consist of the reverse complement of the inline barcode(s) followed by the read 2 adapter
    """
    chunks = []
    for barcode_id in barcode_ids:
        full = revcomp(PE_BARCODES[barcode_id]) + R2_ADAPTER
        for chunk in sliding_windows(full):
            if chunk not in chunks:
                chunks.append(chunk)

    return chunks


def pe_read1_5p_adapters(barcode_ids):
    """
    Returns 5' adapters of read 1 for paired-end eCLIP (cutadapt -g)
    """
    adapters = []
    for barcode_id in barcode_ids:
        if PE_BARCODES[barcode_id]:
            adapters.append(PE_5P_ADAPTER + PE_BARCODES[barcode_id])

    return adapters


def barcodes_fasta():
    """Returns barcode FASTA (as string) used by eclipdemux (barcode sequence is 5' to 3' of read 1)"""
    lines = []
    for barcode_id, seq in PE_BARCODES.items():
        if seq:
            lines.extend([f">{barcode_id}", seq])

    return "\n".join(lines) + "\n"
