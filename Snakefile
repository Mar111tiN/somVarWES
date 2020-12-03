import yaml

# ############ SETUP ##############################
configfile: "configs/config_devel.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']
# extract the scriptdir for creating shell_paths
snakedir = os.path.dirname(workflow.snakefile)
scriptdir = os.path.join(snakedir, "scripts")


# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"


# retrieve the file_df with all the file paths from the samplesheet
sample_df, short_sample_df = get_files(config['inputdirs'], config['samples']['samplesheet'])
chrom_list = get_chrom_list(config)
TN_list = get_tumor_normal_pairs(sample_df)

# ############ INCLUDES ##############################  
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
    tumor_normal = "[^_/.]+",
    type = "[^_/.]+",
    tumor = "[A-Za-z123]+",
    normal = "[A-Za-z123]+",
    chrom = "(chr)?[0-9XY]+",
    filter = "filter[0-9]+",


# ############## MASTER RULE ##############################################

rule all:
    input:
        expand("filter/{tumor_normal_pair}.filter2.loose.csv", tumor_normal_pair=TN_list),
        expand("filterbam/{tumor_normal_pair}.filter2.IGVnav.txt", tumor_normal_pair=TN_list)


###########################################################################


# print out of installed tools
onstart:
    print("    SOMVARWES PIPELINE STARTING.......")
    print('Variant calling using bam files')
    print('bam:', short_sample_df.loc[:, ['bam_path']])
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("SOMVARWES workflow finished - everything ran smoothly")
    if config['cleanup']:
        split_table_pattern = 'chr[^.]+\.csv'
        shell("rm -rf pileup varscan")
        # remove split table/..chr.csv
        shell("ls table/ egrep '{}' | sed 's_^_table/_' | xargs -r rm -f")

