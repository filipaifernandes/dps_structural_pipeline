FROM snakemake/snakemake:latest

RUN conda install -y -c bioconda -c conda-forge \
    biopython \
    fasttree \
    requests \
    pyyaml \
    python=3.13 \
    pandas \
    matplotlib \
    seaborn \
    && conda clean -afy

WORKDIR /workflow
