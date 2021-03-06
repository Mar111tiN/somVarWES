from yaml import CLoader as Loader, load, dump
from subprocess import run
# ############ SETUP ##############################
configfile: "configs/active_config.yaml"
# configfile: "configs/config.json"

workdir: config['workdir']
# extract the scriptdir for creating shell_paths
snakedir = os.path.dirname(workflow.snakefile)
scriptdir = os.path.join(snakedir, "scripts")


# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"


# load the sample independent config file
config = add_config(config, config_name="general")
# load the CNV config
config = add_config(config, config_name="CNV")


# retrieve the file_df with all the file paths from the samplesheet
sample_df, short_sample_df = get_files(config['inputdirs'], config['samples']['samplesheet'])
chrom_list = get_chrom_list(config)
TN_list = get_tumor_normal_pairs(sample_df)
print(short_sample_df)
print(get_ASCAT_list(TN_list))
# ############ INCLUDES ##############################  
include: "includes/varscan.snk"
include: "includes/annotate.snk"
include: "includes/EB.snk"
include: "includes/HDR.snk"
include: "includes/filter.snk"
include: "includes/CNV.snk"

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
        expand("filterbam/{tumor_normal_pair}.filter2.IGVnav.txt", tumor_normal_pair=TN_list),
        # expand("CNV/{s_tn}.cnv.snp.gz", s_tn=get_STN(TN_list)),
        expand("CNV/{s_tn_t}.tumour.png", s_tn_t=get_ASCAT_list(TN_list)),
        expand("SNPlot/{tumor_normal_pair}.SNPmatch.jpg", tumor_normal_pair=TN_list)

###########################################################################


# print out of installed tools
onstart:
    show_output("    SOMVARWES PIPELINE STARTING.......", time=True)

    print('bam:', short_sample_df.loc[:, ['bam_path']])
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # # touch all conda env files in the workdir so the scratch-cleanup does not clean them (BIH-cluster specific workaround)
    # run(f"find {config['workdir']}/.snakemake -type f -exec touch {{}} +", shell=True)
    show_output("SOMVARWES workflow finished - everything ran smoothly", time=True, color="success")
    # if config['cleanup']:
    #     split_table_pattern = 'chr[^.]+\.csv'
    #     shell("rm -rf pileup varscan")
    #     # remove split table/..chr.csv
    #     shell("ls table/ egrep '{}' | sed 's_^_table/_' | xargs -r rm -f")

