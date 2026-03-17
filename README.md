# dps_structural_pipeline
export MODELLER_KEY=YOUR_MODELLER_LICENSE_HERE
snakemake --use-singularity --cores 4


export MODELLER_KEY=YOUR_MODELLER_LICENSE_HERE
export APPTAINERENV_MODELLER_KEY=$MODELLER_KEY
snakemake --use-singularity --cores 4
