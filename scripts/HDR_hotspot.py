import os
import pandas as pd
from multiprocessing import Pool
import numpy as np
import math
from io import StringIO
from functools import partial
from subprocess import Popen, PIPE, run
from script_utils import show_command, show_output
from time import sleep

# ############ FILTER_BAM UTILS ########################################################


def reduce_regions(df, padding):
    '''
    takes a mutation list and returns a region list using padding
    overlapping regions are reduced to one using the gap strategy
    '''

    df = df.sort_values('Start')
    df['Start'] = df['Start'] - padding
    df['End'] = df['End'] + padding
    # find the break points
    # if Start is greater than previous End (using shift), this is a gap --> df['gap'] = 1
    df['gap'] = df['Start'].gt(df['End'].shift()).astype('int')
    # id different reads according to gap
    # cumulative sum does not increase at df['gap'] == 0 and so these consecutive stretches are grouped together
    df['gap'] = df['gap'].cumsum()
    # groupby the coverage break group and condense individual coverage islands
    # agg has to contain the neccessary shared columns TransLength because it is needed for coverage computation
    df = df.groupby('gap').agg({'Chr': 'first', 'Start': 'min', 'End': 'max'})
    return df.reset_index('gap').drop(columns='gap')


def mut2bed(mut_df, chrom, padding, output):
    # check for filter_bam folder (using general declaration for easy folder name changing)
    folder = os.path.split(output)[0]
    if not os.path.isdir(folder):
        os.makedirs(folder)
    # empty file
    if not len(mut_df.index):
        mut_df.to_csv(output, index=False, sep='\t', header=False)
        return output
    # get the bedfile with padded and collapsed regions
    bed_df = reduce_regions(mut_df.sort_values(
        ['Chr', 'Start']).iloc[:, :5], padding)

    # write bed_df to file
    bed_df.to_csv(output, index=False, sep='\t', header=False)
    return output


def bam2hotspot(bam_file, chrom, HDR_config, mut_df):
    chrom_seq = os.path.join(HDR_config['genome_split_path'], f"{chrom}.fa")

    # unwrap the mawk tools
    min_Alt = HDR_config['minAltSum']
    pile2hotspot = HDR_config['pile2hotspot']

    # assign bedfile name
    proc = f".{os.getpid()}" if HDR_config['multi'] else ""
    bed_file = bam_file.replace(".bam", f".{chrom}{proc}.bed")

    # create the bedfile for the mpileup command
    bed_file = mut2bed(
        mut_df, chrom, padding=HDR_config['PAD'], output=bed_file)
    print(bed_file, os.path.isfile(bed_file))
    sleep(2)
    pileup_cmd = f"samtools mpileup -l {bed_file} -f {chrom_seq} -q {HDR_config['MINq']} -Q {HDR_config['MINQ']} {bam_file}"

    cmd = f"{pileup_cmd} | {pile2hotspot} {min_Alt}"
    show_command(cmd)
    hotspot_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    run(f"rm {bed_file}", shell=True)
    return hotspot_df


def bam2hotspot_multi(bam_file, chrom, HDR_config, mut_df, threads):
    # compute true HDRs using multicores
    # MULTITHREADING
    # use half the threads because 2 processes are required for bam2hotspot
    threads = int(threads/2)
    show_output(f"Using {threads} cores for hotspot detection on {chrom}")
    hot_pool = Pool(threads)
    # split
    # minimal length of 10 mutations per thread
    MINLEN = 10
    split_factor = min(math.ceil(len(mut_df.index) / MINLEN), threads)
    mut_split = np.array_split(mut_df, split_factor)
    # apply
    HDR_config['multi'] = True
    hot_dfs = hot_pool.map(
        partial(bam2hotspot, bam_file, chrom, HDR_config), mut_split)
    hot_pool.close()
    hot_pool.terminate()
    hot_pool.join()
    # concat and return
    # there is a chance for overlapping
    hot_df = pd.concat(hot_dfs).drop_duplicates()
    return hot_df


def pileup2hotspot(mut_df, pileup_file, chrom, HDR_config):
    # unwrap the mawk tool
    pile2hotspot = HDR_config['pile2hotspot_chrom']
    min_Alt = HDR_config['minAltSum']
    cmd = f"cat {pileup_file} | {pile2hotspot} {chrom} {min_Alt}"
    show_command(cmd)
    hotspot_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return hotspot_df


# LEGACY ###################################################
def get_count_pileup(filename, sample_count=1):
    '''
    get_count_pileup for different numbers of pileups
    '''
    cols = [0, 1, 2]
    cols += [i + 3 * sample_count for i in cols]
    pileup = pd.read_csv(
        filename,
        sep='\t',
        header=None,
        usecols=cols,
        names=['Chr', 'Pos', 'Ref', 'Depth', 'read', 'Qual']
    )

    pileup['read'] = pileup['read'].str.upper()
    for base in 'ACTG':
        pileup[base] = pileup['read'].str.count(base)
    pileup['AltSum'] = pileup[['A', 'C', 'T', 'G']].max(axis=1)
    pileup['AltRatio'] = pileup['AltSum'] / pileup['Depth']
    return pileup


def filter_hotspots(pileup_df, HDR_config):
    '''
    filters from the pileup_file the putative hotspots of HDR
    '''

    minAlt = HDR_config['minAltSum']
    minRatio = HDR_config['minAltRatio']
    maxRatio = HDR_config['maxAltRatio']
    hotspot_df = pileup_df.query(
        '(AltSum >= @minAlt) and (@minRatio <= AltRatio <= @maxRatio)')
    return hotspot_df


def get_hotspot_df(pileup_file, HDR_config):
    '''
    here comes the pileup computation straight from the bam file
    '''
    pileup_df = get_count_pileup(pileup_file)
    show_output(f"Loading pileup file {pileup_file} finished.")
    hotspot_df = filter_hotspots(pileup_df, HDR_config)
    return hotspot_df
