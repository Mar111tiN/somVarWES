#!/bin/bash

# can be overruled on CLI with -J <NAME>
#SBATCH --job-name DEVEL
#SBATCH --output=logs/%x-%j.log
#SBATCH -n=2
#SBATCH --nodes=1   
#SBATCH --partition=medium
#SBATCH --time=10:00:00


export LOGDIR=logs/${SLURM_JOB_NAME}/${SLURM_JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp
# export WRKDIR=$HOME/work/projects/whWES

mkdir -p $LOGDIR
unset DRMAA_LIBRARY_PATH

# somehow my environments are not set
# have to set it explicitly
# conda activate somVar-EB
conda activate WES-env
# outputs every output to the terminal
set -x

# !!! leading white space is important
DRMAA=" -p medium -t 10:00:00 --mem-per-cpu=3500M --nodes=1 -n {threads}"
DRMAA="$DRMAA -o $LOGDIR/%x-%j.log"
snakemake --unlock --rerun-incomplete
snakemake --dag | dot -Tsvg > dax/dag.svg
snakemake --use-conda  --rerun-incomplete --restart-times 3 --drmaa "$DRMAA" -prkj 4
# -k ..keep going if job fails
# -p ..print out shell commands
