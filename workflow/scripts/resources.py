import os


class Resources:
    """Gets URLs and file names of the reference files used by the ENCODE eCLIP pipeline
    for a given genome.

    ENCODE eCLIP uses UCSC-style chromosome names (chr1, chr2, ...) throughout, so
    genome and annotation are from GENCODE (not Ensembl). The GENCODE releases below
    match the annotations that are built in to CLIPper (see `clipper_species`).
    """

    # create genome directory
    os.makedirs("resources/", exist_ok=True)

    def __init__(self, genome):
        self.genome = genome

        if genome == "hg38":
            base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/"
            self.fasta_url = f"{base_url}GRCh38.primary_assembly.genome.fa.gz"
            self.gtf_url = f"{base_url}gencode.v29.annotation.gtf.gz"
            self.clipper_species = "GRCh38_v29e"

            # ENCODE eCLIP blacklist (GRCh38)
            self.blacklist_url = "https://www.encodeproject.org/files/ENCFF269URO/@@download/ENCFF269URO.bed.gz"
            self.blacklist_stranded = True  # BED6

            # Dfam clade and rDNA repeating unit (NCBI) for the repeat element reference
            self.dfam_clade = "Homo sapiens"
            self.rdna_accession = "U13369.1"

        elif genome == "hg19":
            base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/"
            self.fasta_url = f"{base_url}GRCh37.p13.genome.fa.gz"
            self.gtf_url = f"{base_url}gencode.v19.annotation.gtf.gz"
            self.clipper_species = "hg19"

            # ENCODE eCLIP blacklist (hg19)
            self.blacklist_url = "https://www.encodeproject.org/files/ENCFF039QTN/@@download/ENCFF039QTN.bed.gz"
            self.blacklist_stranded = True  # BED6

            self.dfam_clade = "Homo sapiens"
            self.rdna_accession = "U13369.1"

        elif genome == "mm10":
            base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"
            self.fasta_url = f"{base_url}GRCm38.primary_assembly.genome.fa.gz"
            self.gtf_url = f"{base_url}gencode.vM25.annotation.gtf.gz"
            self.clipper_species = "mm10v25"

            # ENCODE exclusion list (mm10)
            self.blacklist_url = "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
            self.blacklist_stranded = False  # BED3

            self.dfam_clade = "Mus musculus"
            self.rdna_accession = "BK000964.3"

        else:
            raise ValueError(f"Genome {genome} not supported (hg38, hg19 or mm10)")

        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        self.blacklist = self._file_from_url(self.blacklist_url)

        # generated files
        self.fai = f"{self.fasta}.fai"
        self.chrom_sizes = f"resources/{genome}_chrom.sizes"
        self.repeat_fasta = f"resources/{genome}_repeat_elements.fa"
        self.star_index = f"resources/star_index/{genome}"
        self.star_repeat_index = f"resources/star_index/{genome}_repeat_elements"

    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file"""
        return f"resources/{os.path.basename(url).replace('.gz','')}"
