#!/bin/bash

# can be overruled on CLI with -J <NAME>
#SBATCH --job-name=somVar

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=slogs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=2

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=500M

# Set the partition to be used (-p or --partition).
#SBATCH --partition=long

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=20:00:00

SNAKE_HOME=$(pwd);

export LOGDIR=${SNAKE_HOME}/slogs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp;
mkdir -p $LOGDIR;

set -x;

unset DRMAA_LIBRARY_PATH

# make conda available
eval "$($(which conda) shell.bash hook)"
# activate snakemake env
conda activate WES-env;
echo $CONDA_PREFIX "activated";


# !!! leading white space is important
DRMAA=" -p {cluster.partition} -t {cluster.t} --mem-per-cpu={cluster.mem} -J {cluster.name} --nodes={cluster.nodes} -n {cluster.threads}";
DRMAA="$DRMAA -o ${LOGDIR}/{rule}-%j.log";
snakemake --unlock -s SnakefileRerun --rerun-incomplete;
snakemake --dag | awk '$0 ~ "digraph" {p=1} p' | dot -Tsvg > dax/dag.svg;
snakemake --use-conda --rerun-incomplete --cluster-config configs/cluster/somvar-cluster.json --drmaa "$DRMAA" -prk -j 1000;
# -k ..keep going if job fails
# -p ..print out shell commands
