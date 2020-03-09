import pandas as pd
import os
from HDR_utils import masterHDR
from script_utils import show_output


# ############## SNAKEMAKE ################################
w = snakemake.wildcards
config = snakemake.config
params = snakemake.params
MIN_SIM = params.min_sim
PAD = min(config['HDR']['padding'], config['filter_bam']['padding'])


i = snakemake.input


tumor_bam = os.path.split(i.filter_bam)[1]
# remove the normal part of the tumor-normal descriptor
tumor_bam = tumor_bam.replace(
    f"-{config['samples']['normal'][0]}", '').replace('.done', '.bam')
# folder comes from config
tumor_bam = os.path.join(config['filter_bam']['folder'], tumor_bam)


filter_file = i.filter_file
filter_pileup = i.pileup
out_file = snakemake.output.HDR


show_output(
    f'Starting HDR analysis of {filter_file}. [MIN_SIM={MIN_SIM}, PAD={PAD}]')
# GET THE mutation file for masterHDR
show_output(f'Importing {filter_file} for HDR detection', time=False)
filter_df = pd.read_csv(filter_file, sep='\t').loc[:, [
    'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']]

HDR_df = masterHDR(
    pileup_file=filter_pileup,
    bam_file=tumor_bam,
    filter_df=filter_df,
    min_sim=MIN_SIM,
    padding=PAD
)

HDR_df['HDRmeanSimilarity'] = HDR_df['HDRmeanSimilarity'].round(2)

HDR_len = len(HDR_df.query('HDRcount > 0').index)

show_output(
    f"Found {HDR_len} possible HDR mutations. Writing to file {out_file}")
HDR_df.to_csv(out_file, sep='\t', index=False)
