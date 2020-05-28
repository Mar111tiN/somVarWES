#!/bin/bash

#$ -V
#$ -j y
#$ -N FCAHsub
#$ -o logs
#$ -r yes
#$ -cwd
#$ -S /bin/bash
#$ -P control
#$ -pe smp 2
#s -l h_rt=10:00:00

export LOGDIR=/fast/users/${USER}/scratch/lox/WES/${JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp
# export WRKDIR=$HOME/work/projects/whWES

mkdir -p $LOGDIR

# somehow my environments are not set
# have to set it explicitly
# conda activate somVar-EB
conda activate WES-env
# outputs every output to the terminal
set -x

# !!! leading white space is important
DRMAA=" -P medium -pe smp {threads}  -l h_rt=01:30:00 -l h_vmem=3.5g"
DRMAA="$DRMAA -V -o $LOGDIR/ -j yes"
snakemake --snakefile Snakefile1 --unlock --rerun-incomplete
snakemake --snakefile Snakefile1 --dag | dot -Tsvg > dax/dag.svg
snakemake --snakefile Snakefile1 --use-conda  --rerun-incomplete --drmaa "$DRMAA" -j 1000 -p -r -k
# -k ..keep going if job fails
# -p ..print out shell commands