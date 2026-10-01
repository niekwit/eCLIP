FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2b918a2150870c895ad2fe70382897809d83c86abde4d06f7d2a62090b9ae2be"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/7cdb9e3180904a671be6b170bf7651d2
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.7.2
RUN mkdir -p /conda-envs/7cdb9e3180904a671be6b170bf7651d2
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/fastqc/environment.yaml /conda-envs/7cdb9e3180904a671be6b170bf7651d2/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/84188cb48cfebe3cadc377f2f5001295
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.22
#     - snakemake-wrapper-utils =0.7.2
RUN mkdir -p /conda-envs/84188cb48cfebe3cadc377f2f5001295
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/faidx/environment.yaml /conda-envs/84188cb48cfebe3cadc377f2f5001295/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/93bab9ed6097f300f8667bc95cf9fd61
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.22
RUN mkdir -p /conda-envs/93bab9ed6097f300f8667bc95cf9fd61
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/index/environment.yaml /conda-envs/93bab9ed6097f300f8667bc95cf9fd61/environment.yaml

# Conda environment:
#   source: workflow/envs/bigwig.yaml
#   prefix: /conda-envs/e2b1fb9dc0089fc6acc826cfb3c22d05
#   name: bigwig
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.12
#     - pysam=0.24.1
#     - bedtools=2.31.1
#     - samtools=1.21
#     - ucsc-bedgraphtobigwig=482
RUN mkdir -p /conda-envs/e2b1fb9dc0089fc6acc826cfb3c22d05
COPY workflow/envs/bigwig.yaml /conda-envs/e2b1fb9dc0089fc6acc826cfb3c22d05/environment.yaml

# Conda environment:
#   source: workflow/envs/clipper.yaml
#   prefix: /conda-envs/addedf23fe3ece2b9c64f11c291df800
#   name: clipper
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     # CLIPper itself is installed by clipper.post-deploy.sh
#     - python=3.7
#     - pip
#     - git
#     - cxx-compiler
#     - cython=0.29
#     - numpy=1.21
#     - scipy=1.5
#     - pandas=1.0
#     - scikit-learn=0.23
#     - matplotlib-base=3.2
#     - pysam=0.15
#     - pybedtools=0.8.1
#     - htseq=0.11
#     - bedtools=2.29
#     - samtools=1.9
RUN mkdir -p /conda-envs/addedf23fe3ece2b9c64f11c291df800
COPY workflow/envs/clipper.yaml /conda-envs/addedf23fe3ece2b9c64f11c291df800/environment.yaml

# Conda environment:
#   source: workflow/envs/idr.yaml
#   prefix: /conda-envs/8139c61937c457a62159880ff4124b7a
#   name: idr
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     # ENCODE used IDR 2.0.2 (not installable anymore, requires Python 3.5)
#     - idr=2.0.4.2
#     - python=3.10
#     - numpy<1.24
RUN mkdir -p /conda-envs/8139c61937c457a62159880ff4124b7a
COPY workflow/envs/idr.yaml /conda-envs/8139c61937c457a62159880ff4124b7a/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yaml
#   prefix: /conda-envs/12f97f22c7e34ec074626c2614a31e37
#   name: mapping
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - star=2.7.11b
#     - samtools=1.21
#     - fastq-tools=0.8.3
#     - pigz=2.8
#     - gawk=5.3.1
#     - wget=1.21.4
#     - multiqc=1.25.2
#     - python=3.12
RUN mkdir -p /conda-envs/12f97f22c7e34ec074626c2614a31e37
COPY workflow/envs/mapping.yaml /conda-envs/12f97f22c7e34ec074626c2614a31e37/environment.yaml

# Conda environment:
#   source: workflow/envs/peaks.yaml
#   prefix: /conda-envs/d84036563bf2918d3c24ab0b0eb188e9
#   name: peaks
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - perl=5.32
#     - perl-statistics-basic=1.6611
#     - perl-statistics-distributions=1.02
#     - perl-statistics-r=0.34
#     - r-base=4.3
#     - python=3.10
#     - pandas=2.1.4
#     - samtools=1.21
#     - bedtools=2.31.1
#     - ucsc-bedtobigbed=482
RUN mkdir -p /conda-envs/d84036563bf2918d3c24ab0b0eb188e9
COPY workflow/envs/peaks.yaml /conda-envs/d84036563bf2918d3c24ab0b0eb188e9/environment.yaml

# Conda environment:
#   source: workflow/envs/te_repeats.yaml
#   prefix: /conda-envs/c56cc2c1aa93a54ea1d2306495a932af
#   name: te_repeats
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.12
#     - pandas=2.2.3
#     - scipy=1.14.1
#     - statsmodels=0.14.4
#     - bedtools=2.31.1
RUN mkdir -p /conda-envs/c56cc2c1aa93a54ea1d2306495a932af
COPY workflow/envs/te_repeats.yaml /conda-envs/c56cc2c1aa93a54ea1d2306495a932af/environment.yaml

# Conda environment:
#   source: workflow/envs/trimming.yaml
#   prefix: /conda-envs/8cdf60bd329f5889e9ddb1100adf2930
#   name: trimming
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - cutadapt=5.1
#     - pigz=2.8
RUN mkdir -p /conda-envs/8cdf60bd329f5889e9ddb1100adf2930
COPY workflow/envs/trimming.yaml /conda-envs/8cdf60bd329f5889e9ddb1100adf2930/environment.yaml

# Conda environment:
#   source: workflow/envs/yeolab_py2.yaml
#   prefix: /conda-envs/f0968f5032e34a9492efe3651d397aff
#   name: yeolab_py2
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     # eclipdemux and barcodecollapsepe.py (Yeo lab) are Python 2 code (used with Python 2.7 by ENCODE)
#     - python=2.7
#     - pip
#     - pysam=0.15.4
#     - future=0.18.2
#     - samtools=1.9
#     - pip:
#       - git+https://github.com/YeoLab/eclipdemux.git@2ccc6b6b15f92b4613bb489c693d0567313baf3b
RUN mkdir -p /conda-envs/f0968f5032e34a9492efe3651d397aff
COPY workflow/envs/yeolab_py2.yaml /conda-envs/f0968f5032e34a9492efe3651d397aff/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/7cdb9e3180904a671be6b170bf7651d2 --file /conda-envs/7cdb9e3180904a671be6b170bf7651d2/environment.yaml && \
    conda env create --prefix /conda-envs/84188cb48cfebe3cadc377f2f5001295 --file /conda-envs/84188cb48cfebe3cadc377f2f5001295/environment.yaml && \
    conda env create --prefix /conda-envs/93bab9ed6097f300f8667bc95cf9fd61 --file /conda-envs/93bab9ed6097f300f8667bc95cf9fd61/environment.yaml && \
    conda env create --prefix /conda-envs/e2b1fb9dc0089fc6acc826cfb3c22d05 --file /conda-envs/e2b1fb9dc0089fc6acc826cfb3c22d05/environment.yaml && \
    conda env create --prefix /conda-envs/addedf23fe3ece2b9c64f11c291df800 --file /conda-envs/addedf23fe3ece2b9c64f11c291df800/environment.yaml && \
    conda env create --prefix /conda-envs/8139c61937c457a62159880ff4124b7a --file /conda-envs/8139c61937c457a62159880ff4124b7a/environment.yaml && \
    conda env create --prefix /conda-envs/12f97f22c7e34ec074626c2614a31e37 --file /conda-envs/12f97f22c7e34ec074626c2614a31e37/environment.yaml && \
    conda env create --prefix /conda-envs/d84036563bf2918d3c24ab0b0eb188e9 --file /conda-envs/d84036563bf2918d3c24ab0b0eb188e9/environment.yaml && \
    conda env create --prefix /conda-envs/c56cc2c1aa93a54ea1d2306495a932af --file /conda-envs/c56cc2c1aa93a54ea1d2306495a932af/environment.yaml && \
    conda env create --prefix /conda-envs/8cdf60bd329f5889e9ddb1100adf2930 --file /conda-envs/8cdf60bd329f5889e9ddb1100adf2930/environment.yaml && \
    conda env create --prefix /conda-envs/f0968f5032e34a9492efe3651d397aff --file /conda-envs/f0968f5032e34a9492efe3651d397aff/environment.yaml && \
    conda clean --all -y

# Step 4: Run post-deploy scripts (not added by snakemake --containerize)

# CLIPper is not on bioconda: build it from source into the clipper environment
COPY workflow/envs/clipper.post-deploy.sh /tmp/clipper.post-deploy.sh
RUN conda run -p /conda-envs/addedf23fe3ece2b9c64f11c291df800 bash /tmp/clipper.post-deploy.sh && \
    /conda-envs/addedf23fe3ece2b9c64f11c291df800/bin/clipper --help > /dev/null && \
    rm /tmp/clipper.post-deploy.sh
