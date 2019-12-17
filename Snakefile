import os
import re
import yaml
import argparse
import math

# ############ SETUP ##############################
configfile: "configs/config_devel.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']

# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"

# FILE SEARCH ###########
#   included_files = tumor_normal_pairs = sample_types = None
#  if not included_files:  # do this only once
#  included_files, tumor_normal_pairs, sample_types = get_files(config['inputdir'])
# included files : list of (fastq_path, sample, type, read, tumor/normal) tuples
# tumor_normal_pairs : list of 'sample_tumor_normal' strings
# sample_types : list of 'sample_type' strings


# retrieve the file_df with all the file paths from the samplesheet
sample_df = get_files(config['inputdirs'], config['samples']['samplesheet'])
chrom_list = get_chrom_list(config)


# ############ INCLUDES ##############################
include: "includes/fastq.snk"
include: "includes/QC.snk"
include: "includes/ubam.snk"
include: "includes/map.snk"
include: "includes/splitBAM.snk"
include: "includes/processBAM.snk"
include: "includes/dedup.snk"
include: "includes/varscan.snk"
include: "includes/annotate.snk"
include: "includes/EB.snk"
include: "includes/filter.snk"

# convenience variables
ref_gen = full_path('genome')
# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^_/.]+",
    type = "[^_/.]+",
    read = "[^_/.]+",
    tumor_norm = "[^_/.]+",
    tumor = "[A-Za-z]+",
    norm = "[A-Za-z]+",
    split = "[0-9]+",
    read_or_index = "[^_/.]+",
    trim = "[^_/.]+",
    chrom = "[^_/.]*[0-9XY]+",
    filter = "[A-Za-z]+"


# extract the filter list for active filters
active_filter_list = [f for f in config['filter']['filters'].keys() if config['filter']['filters'][f]['run']]


# ############## MASTER RULE ##############################################

rule all:
    input:
        "QC/fastQC.html",
        "QC/libraryQC.html",
        "QC/insertQC.html",
        expand("coverBED/{samples}.txt", samples=sample_df.index)
        # expand("filter/{tumor_normal_pair}.{filter}.csv", tumor_normal_pair=get_tumor_normal_pairs(sample_df), filter=active_filter_list),
        # expand("filter_bam/filter_bam/{tumor_normal_pair}.{filter}.bam", tumor_normal_pair=get_tumor_normal_pairs(sample_df), filter=active_filter_list)

###########################################################################


# print out of installed tools
onstart:
    print("    EXOM SEQUENCING PIPELINE STARTING.......")
    print('samples', sample_df['R1'])
    ##########################
    # shell("echo Conda-environment: $CONDA_PREFIX")
    # shell('echo $PATH')
    # write config to the results directory
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder
#     shell("conda list | show_awk")
#     shell("ls -l ${{TOOLS}}/bulltools/links | sed -nr '/->/s_.* ([^/]+ -> .*)_  \\1_p' ")
    # create scratch/tmp/hWES folder for storing stuff


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")
