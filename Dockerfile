FROM snakemake/snakemake:7.32.4

RUN conda install -y -c bioconda -c conda-forge \
    biopython \
    fasttree \
    requests \
    pyyaml \
    python=3.13 && \
    conda clean -afy

WORKDIR /workflow
