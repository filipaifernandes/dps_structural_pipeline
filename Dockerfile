FROM snakemake/snakemake:latest

RUN conda install -c bioconda -c conda-forge \
    biopython \
    fasttree \
    requests \
    python=3.13

RUN conda install -c salilab -c conda-forge modeller

WORKDIR /workflow
