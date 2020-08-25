import pandas as pd
from HDR_core import get_HDR_multi, get_filter_hdr_multi
from script_utils import show_output
from HDR_hotspot import bam2hotspot_multi, pileup2hotspot


COLS = [
    'Chr',
    'Start',
    'End',
    'Ref',
    'Alt',
    'Gene',
    'HDRcand',
    'HDRcount',
    'HDRinfo',
]


def HDR_master(mut_file, bam_file, chrom, threads, HDR_config, pileup_file=''):
    '''
    collects the imput and 
    '''

    HDR_config['multi'] = (threads > 1)
    # Get the mutation file for masterHDR
    show_output(f'Importing {mut_file} for HDR detection', time=False)
    mut_df = pd.read_csv(mut_file, sep='\t').loc[:, [
        'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']]
    # check for empty file
    if len(mut_df.index) == 0:
        # return an empty df
        show_output(
            f"Mutation file {mut_file} is empty! Writing empty file", color="warning")
        return pd.DataFrame(columns=COLS)

    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ['chrX', 'chrY']
    mut_df['Chr'] = pd.Categorical(mut_df['Chr'], chrom_list)
    mut_df = mut_df.sort_values(['Chr', 'Start'])

    if chrom:
        HDR_df = HDR_run(mut_df, mut_file, bam_file, chrom,
                         threads, HDR_config, pileup_file)
    else:
        # if no chroms are given
        chroms = mut_df['Chr'].unique()
        HDR_dfs = []
        for chrom in chroms:
            hdr_df = HDR_run(mut_df, mut_file, bam_file, chrom,
                             threads, HDR_config, pileup_file)
            HDR_dfs.append(hdr_df)

        HDR_df = pd.concat(HDR_dfs).sort_values('Chr')

    # merge the HDR output into mut_df
    HDR_df = mut_df.merge(
        HDR_df, on=COLS[:6], how="outer").sort_values(COLS[:2])
    for col in COLS[6:8]:
        HDR_df[col] = HDR_df[col].fillna(0).astype(int)
    HDR_df['HDRinfo'] = HDR_df['HDRinfo'].fillna('no HDR in vincinity')

    return HDR_df


def HDR_run(mut_df, mut_file, bam_file, chrom, threads, HDR_config, pileup_file):
    '''
    does all the work and gets directly called from snakemake/CLI wrapper
    '''
    if HDR_config['multi']:
        cores = f" using {threads} cores"
    show_output(
        f"Starting HDR analysis of {mut_file} on {chrom}{cores}. [MIN_SIM={HDR_config['MINSIM']}, PAD={HDR_config['PAD']}]")
    # restrict mut_df to chrom
    mut_df = mut_df.query('Chr == @chrom')

    ############ FIND HOTSPOTS ###########################
    if pileup_file:
        source = pileup_file
        hotspot_df = pileup2hotspot(
            mut_df, pileup_file, chrom=chrom, HDR_config=HDR_config)
    else:
        source = bam_file
        # hotspot_df = bam2hotspot(bam_file, chrom=chrom, HDR_config=HDR_config, mut_df=mut_df)
        hotspot_df = bam2hotspot_multi(
            bam_file, chrom=chrom, HDR_config=HDR_config, mut_df=mut_df, threads=threads)
    hotspot_len = len(hotspot_df.index)

    # check empty
    if hotspot_len == 0:
        show_output(
            f"Found no putative HDR lanes in {source}.")
        return pd.DataFrame(columns=COLS)

    show_output(
        f"Detected {hotspot_len} putative HDR lanes in {source}.")

    ############ FILTER HOTSPOTS ###########################
    filter_HDR = get_filter_hdr_multi(mut_df, hotspot_df, threads, HDR_config)
    filter_HDR_len = len(filter_HDR.index)

    if filter_HDR_len == 0:
        show_output(
            f"Found no HDR-rich mutations in {source}", multi=False)
        return pd.DataFrame(columns=COLS)
    show_output(
        f"Found a total of {filter_HDR_len} HDR-rich mutations", multi=False)

    HDR_df = get_HDR_multi(filter_HDR, bam_file=bam_file, hotspot_df=hotspot_df, threads=threads,
                           chrom=chrom, HDR_config=HDR_config)

    return HDR_df
