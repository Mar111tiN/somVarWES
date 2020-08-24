import pandas as pd
from script_utils import show_output, show_command
from multiprocessing import Pool
import numpy as np
import math
from functools import partial
from io import StringIO
from subprocess import PIPE, run
from HDR_bam import bam2df2
from HDR_info import get_HDR_info, get_HDR_count


def get_filter_hdr(hotspot_df, HDR_config, hdr_df):
    # enumerate the HDRs in vicinity of mutations
    hdr_df.loc[:, 'HDRcand'] = hdr_df.apply(
        get_HDR_count, axis=1, args=(hotspot_df,), padding=HDR_config['PAD'])
    min_HDR_count = HDR_config['MinHDRCount']
    # filter HDR-rich mutations
    filter_HDR = hdr_df.query('HDRcand >= @min_HDR_count')
    show_output(
        f"Found {len(filter_HDR.index)} HDR-rich mutations", multi=True)
    return filter_HDR


def get_filter_hdr_multi(HDR_df, hotspot_df, threads, HDR_config):
    # compute the filter_hdr using multicores
    # MULTITHREADING
    # pool
    show_output(f"Using {threads} cores for filtering HDR-rich mutations")
    filter_pool = Pool(threads)
    # split
    # minimal length of 200 lines
    split_factor = min(math.ceil(len(HDR_df.index) / 200), threads)
    hdr_split = np.array_split(HDR_df, split_factor)
    # apply
    filter_HDRs = filter_pool.map(
        partial(get_filter_hdr, hotspot_df, HDR_config), hdr_split)
    # concat
    filter_HDR = pd.concat(filter_HDRs).sort_values(['Chr', 'Start'])
    filter_pool.close()
    filter_pool.terminate()
    filter_pool.join()
    return filter_HDR


def get_HDR(bam_file, hotspot_df, chrom, HDR_config, filter_HDR):

    # get the bam_df for analysis in single-read resolution
    # get the required region for bam analysis
    padding = HDR_config['PAD']
    _min = filter_HDR['Start'].min() - padding
    _max = filter_HDR['Start'].max() + padding
    region = f"{chrom}:{_min}-{_max}"
    bam_df = bam2df2(bam_file, region=region, HDR_config=HDR_config)
    show_output(
        f'Loaded the bam_file {bam_file} on region {region} for read analysis - {len(bam_df.index)}', multi=True)

    filter_HDR[['HDRcount', 'HDRinfo']] = filter_HDR.apply(
        get_HDR_info, axis=1, args=(hotspot_df, bam_df), HDR_config=HDR_config)
    filter_HDR.loc[:, 'HDRcount'] = filter_HDR['HDRcount'].fillna(
        0).astype(int)
    filter_HDR.loc[:, 'HDRinfo'] = filter_HDR['HDRinfo'].fillna('no HDR')
    show_output(
        f"Done analysing {len(filter_HDR.index)} HDR-rich mutations", multi=True)
    return filter_HDR


def get_HDR_multi(filter_HDR, bam_file, hotspot_df, threads,
                  chrom, HDR_config):
    # compute true HDRs using multicores
    # MULTITHREADING
    show_output(f"Using {threads} cores for HDR computation")
    hdr_pool = Pool(threads)
    # split
    # minimal length of 2 mutations
    split_factor = min(math.ceil(len(filter_HDR.index) / 2), threads)
    filter_hdr_split = np.array_split(filter_HDR, split_factor)
    # apply
    HDR_dfs = hdr_pool.map(partial(get_HDR, bam_file, hotspot_df,
                                   chrom, HDR_config), filter_hdr_split)
    hdr_pool.close()
    hdr_pool.terminate()
    hdr_pool.join()
    # concat and return
    return pd.concat(HDR_dfs, sort=True)
