FROM snakemake/snakemake:latest

RUN conda install -y -c bioconda -c conda-forge \
    biopython \
    iqtree \
    requests \
    pyyaml \
    python \
    pandas \
    matplotlib \
    seaborn \
    && conda clean -afy

WORKDIR /workflow
