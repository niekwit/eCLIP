"""
Creates repeat element reference (FASTA) for the repeat element filtering step of ENCODE eCLIP.

ENCODE used RepBase (which is not freely available). This uses the curated consensus sequences
of the species from Dfam (https://www.dfam.org) and the rDNA repeating unit from NCBI as a substitute.
"""

import json
import sys
import time
import urllib.parse
import urllib.request

DFAM_API = "https://www.dfam.org/api/families"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAGE_SIZE = 1000  # maximum allowed by Dfam API


def fetch(url, retries=5):
    """Returns response body of url as string (retries on failure)"""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read().decode()
        except Exception as e:
            print(f"Failed to fetch {url}: {e} (attempt {attempt + 1})", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Could not fetch {url}")


def dfam_url(clade, start, fmt, limit):
    query = urllib.parse.urlencode(
        {
            "format": fmt,
            "clade": clade,
            "clade_relatives": "ancestors",
            "start": start,
            "limit": limit,
        },
        quote_via=urllib.parse.quote,  # Dfam does not accept + for spaces
    )
    return f"{DFAM_API}?{query}"


def dfam_fasta(clade):
    """Returns FASTA of all curated Dfam families of a clade (including families of ancestral clades)"""
    total = json.loads(fetch(dfam_url(clade, 0, "summary", 1)))["total_count"]
    print(f"Fetching {total} Dfam families for {clade}", flush=True)

    records = []
    for start in range(0, total, PAGE_SIZE):
        body = json.loads(fetch(dfam_url(clade, start, "fasta", PAGE_SIZE)))["body"]
        records.append(body.strip())

    fasta = "\n".join(records) + "\n"

    # Header: >DF000000001.4 MIR -> >DF000000001.4|MIR (STAR only uses the first word of a header)
    lines = []
    for line in fasta.split("\n"):
        if line.startswith(">"):
            line = ">" + "|".join(line[1:].split(maxsplit=1)).replace("/", "_")
        lines.append(line)

    n = sum(1 for x in lines if x.startswith(">"))
    if n != total:
        raise RuntimeError(f"Expected {total} Dfam families, but got {n}")

    return "\n".join(lines)


def rdna_fasta(accession):
    """Returns FASTA of rDNA repeating unit from NCBI"""
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
    )
    fasta = fetch(f"{NCBI_EFETCH}?{query}")
    header, *seq = fasta.strip().split("\n")
    if not header.startswith(">"):
        raise RuntimeError(f"Could not fetch {accession} from NCBI")

    return "\n".join([f">rDNA|{accession}"] + seq) + "\n"


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    sys.stdout = log
    sys.stderr = log

    fasta = dfam_fasta(snakemake.params.clade)
    if not fasta.endswith("\n"):
        fasta += "\n"
    fasta += rdna_fasta(snakemake.params.rdna)

    with open(snakemake.output[0], "w") as f:
        f.write(fasta)
