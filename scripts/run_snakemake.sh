#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=4:00:00
#SBATCH --mem=128GB
#SBATCH -o /home/a1018048/slurm/MFM-223_DHT-RNASeq/%x_%j.out
#SBATCH -e /home/a1018048/slurm/MFM-223_DHT-RNASeq/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Cores
CORES=32
if [ -d "/hpcfs" ]; then
	module load arch/arch/haswell
	module load arch/haswell
	module load modulefiles/arch/haswell
	HPC="/hpcfs"
else
    if [ -d "/fast" ]; then
        HPC=/fast
    else
        exit 1
    fi
fi

## Project Root
PROJ=${HPC}/users/a1018048/MFM-2233_DHT-RNASeq

## The environment containing snakemake
micromamba activate snakemake
cd ${PROJ}

## Run snakemake
snakemake \
  --cores ${CORES} \
  --use-conda \
  --wrapper-prefix 'https://raw.githubusercontent.com/snakemake/snakemake-wrappers/'

## Add files to git
bash ${PROJ}/scripts/update_git.sh
