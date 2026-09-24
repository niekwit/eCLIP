# Test data

Small data set for testing the workflow (used by the GitHub Actions workflow in `.github/workflows/main.yml`):

* `reads/`: simulated single-end eCLIP reads (2 IP replicates, 1 size-matched input) from 28 genes of chr21.
* `resources/`: test genome (hg38 chr21 with the real sequence only around these 28 genes, at the original coordinates), the GENCODE v29 annotation of these genes, and the ENCODE eCLIP blacklist. As these files are present, they are not downloaded.
* `config/config.yaml`: CLIPper is restricted to the test genes (`clipper: extra`), otherwise it would process the whole genome annotation.

Run from this directory:

```bash
snakemake -s ../workflow/Snakefile --use-conda --cores 3
```
