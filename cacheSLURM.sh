#!/bin/bash

# can be overruled on CLI with -J <NAME>
#SBATCH --job-name=EBcache

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=slogs/%x-%j.log

#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=500M
#SBATCH --partition=long
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=1-10:00:00

SNAKE_HOME=$(pwd);

export LOGDIR=${HOME}/scratch/slogs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp;
mkdir -p $LOGDIR;

set -x;

# make conda available
eval "$($(which conda) shell.bash hook)"
# have to set it explicitly
conda activate snake-env;
echo $CONDA_PREFIX "activated";


SLURM_CLUSTER="sbatch -p {cluster.partition} -t {cluster.t} --mem-per-cpu={cluster.mem} -J {cluster.name} --nodes={cluster.nodes} -n {cluster.threads}";
SLURM_CLUSTER="$SLURM_CLUSTER -o ${LOGDIR}/{rule}-%j.log";
snakemake --snakefile cacheSnakefile --unlock
snakemake --dag | awk '$0 ~ "digraph" {p=1} p' | dot -Tsvg > dax/cacheDag.svg
snakemake --snakefile cacheSnakefile --use-conda \
--rerun-incomplete \
--cluster-config configs/cluster/PONcache-cluster.json \
--cluster "$SLURM_CLUSTER" \
-prkj 3000
# -k ..keep going if job fails
# -p ..print out shell commands
# -P medium