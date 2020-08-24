import yaml

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
sample_df, short_sample_df = get_files(config['inputdirs'], config['samples']['samplesheet'])
chrom_list = get_chrom_list(config)
TN_list = get_tumor_normal_pairs(sample_df)

# ############ INCLUDES ##############################  
include: "includes/fastq.snk"
include: "includes/QC.snk"
include: "includes/map.snk"
include: "includes/splitBAM.snk"
include: "includes/processBAM.snk"
include: "includes/dedup.snk"
include: "includes/umi_filter.snk"
include: "includes/varscan.snk"
include: "includes/annotate.snk"
include: "includes/EB.snk"
include: "includes/HDR.snk"
include: "includes/filter.snk"

# convenience variables
ref_gen = full_path('genome')
# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^_/.]+",
    type = "[^_/.]+",
    read = "[^_/.]+",
    tumor_normal = "[^_/.]+",
    tumor = "[A-Za-z123]+",
    normal = "[A-Za-z123]+",
    split = "[0-9]+",
    read_or_index = "[^_/.]+",
    trim = "[^_/.]+",
    chrom = "(chr)?[0-9XY]+",
    filter = "filter[0-9]+",
    chrom_split = "[^_/.]+",
    # folder = "^((?!filter).)*$"


# ############## MASTER RULE ##############################################

rule all:
    input:
        QC_output,
        expand("coverBED/{samples}.txt", samples=sample_df.index),
        expand("table/{tumor_normal_pair}.EB.csv", tumor_normal_pair=TN_list),
        expand("filter/{tumor_normal_pair}.filter2.loose.csv", tumor_normal_pair=TN_list),
        expand("filterbam/{tumor_normal_pair}.filter2.IGVnav.txt", tumor_normal_pair=TN_list)

###########################################################################

# print out of installed tools
onstart:
    print("    EXOM SEQUENCING PIPELINE STARTING.......")
    if config['setup']['rerun']:
        print('Rerun variant calling for existing bam files')
        print('bam:', short_sample_df.loc[:, ['bam_path']])
    else:
        print('fastq:', short_sample_df.loc[:, ['R1', 'R2', 'index']])
    ##########################
    # shell("echo Conda-environment: $CONDA_PREFIX")
    # shell('echo $PATH')
    # write config to the results directory
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")

    # cleanup
    if config['setup']['cleanup']['run']:
        no_bams = config['setup']['cleanup']['keep_only_final_bam']

        split_fastq_pattern = '\.[0-9]+\.fastq.gz'
        split_bam_pattern = 'chr[^.]+\..*'
        split_table_pattern = 'chr[^.]+\.csv'

        shell("rm -rf ubam realigned bam_metrics insert_metrics pileup varscan fastqc mapped bam_done")

        # remove split table/..chr.csv
        shell("ls table/ egrep '{}' | sed 's_^_table/_' | xargs -r rm -f")

        # remove split fastqs
        shell("ls fastq | grep -E '{split_fastq_pattern}' | sed 's_^_fastq/_' | xargs -r rm -f")
        
        # remove split recalib bams
        shell("ls recalib | grep -E '{split_bam_pattern}' | sed 's-^-recalib/-' | xargs -r rm -f")

        if no_bams:
            shell("rm -r bam_merge ")
        else:
            # only remove split_bams in bam_merge
            shell("ls bam_merge | grep -E '{split_bam_pattern}' | sed 's-^-bam_merge/-' | xargs -r rm -f")
