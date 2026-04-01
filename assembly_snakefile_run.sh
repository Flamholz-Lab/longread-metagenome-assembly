#!/bin/sh
#SBATCH --cpus-per-task 8
#SBATCH --gpus=1
#SBATCH --partition=hpc_a10_a
#SBATCH --time 24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jmcdonald@rockefeller.edu
###################

echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

source /ru-auth/local/home/jmcdonald/miniconda3/etc/profile.d/conda.sh
conda activate /ru-auth/local/home/jmcdonald/miniconda3/envs/snakemake

snakemake \
    --snakefile assembly_snakefile \
    --configfile ./assembly_snakefile_config.yaml \
    --executor slurm \
    --jobs 10 \
    --default-resources slurm_partition=hpc_a10_a mem_mb=16000 runtime="24h" \
    --use-conda
