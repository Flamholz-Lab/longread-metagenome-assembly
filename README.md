# longread-metagenome-assembly
Assembly of individual microbial genomes from nanopore sequencing of a metagenomic sample

## Inputs
This repository contains a snakemake pipeline that takes a .fastq file as input (e.g. your raw sequencing data) and runs through several assembly and quality checking steps. 

To run, you’ll need three files in your working directory or path:

`/assembly_snakefile.sh` \
`/assembly_snakefile_config.yaml` \
`/assembly_snakefile_run.sh ` 

First edit the config file to reflect your .fastq file and other desired parameters. The actual snakefile (/assembly_snakefile.sh) should not need to be edited. Note that the snakefile utilizes conda environments that I’ve created under my profile in HPC storage. TBD if others in the lab can access those. If not, you will want to make your own conda environments that have the right packages.

To run on the HPC use the command:

`sbatch assembly_snakefile_run.sh` 

The entire pipeline took about 2 hours using 16 threads for a ~5 GB ONT .fastq file containing genomes from 4 strains. Expect time to scale accordingly.

## Outputs
Converged bins -- .fa files are your genomes!: \
`dastool_output/` 

Taxonomic analysis for converged bins -- the most closely related species to your genomes: \
`kraken2_output/`

Quality and completeness analysis of dastool-converged bins: \
busco_output/batch_summary.txt \
`checkm2_output/`

Calculated depth and coverage from polished MetaFlye assembly: \
`depth.txt` \
`coverage.txt`

Bins from different binners: \
`comebin.tsv` \
`metabat2_cut.tsv` \
`semibin.tsv`  

Relative abundance by coverM, uses filtered.fastq against reference genomes: \
`relative_abundance.tsv`
