import pandas as pd
from script_utils import show_output, show_command
from multiprocessing import Pool
import numpy as np
import math
from functools import partial
from io import StringIO
from subprocess import Popen, PIPE, run


def get_count_pileup(filename, sample_count=1):
    '''
    get_count_pileup for different numbers of pileups
    '''
    cols = [0, 1, 2]
    cols += [i + 3 * sample_count for i in [3, 4, 5]]
    pileup = pd.read_csv(
        filename,
        sep='\t',
        header=None,
        usecols=cols,
        names=['Chr', 'Pos', 'Ref', 'depth', 'read', 'Qual']
    )

    pileup['read'] = pileup['read'].str.upper()
    for base in 'ACTG':
        pileup[base] = pileup['read'].str.count(base)
    pileup['AltSum'] = pileup[['A', 'C', 'T', 'G']].max(axis=1)
    pileup['AltRatio'] = pileup['AltSum'] / pileup['depth']
    return pileup


def filter_hotspots(pileup_df, hotspot_config={
    "minAltSum": 4,
    "minAltRatio": 0.1,
    "maxAltRatio": 0.9
}):
    '''
    filters from the pileup_file the putative hotspots of HDR
    '''

    minAlt = hotspot_config['minAltSum']
    minRatio = hotspot_config['minAltRatio']
    maxRatio = hotspot_config['maxAltRatio']
    hotspot_df = pileup_df.query(
        '(AltSum >= @minAlt) and (@minRatio <= AltRatio <= @maxRatio)')
    return hotspot_df


def bam2df2(bam_file, HDR_config, bamtags="", tool='samtools', mut_row=None, region=''):
    '''
    set the region requires 3 threads
    '''

    if region:
        mut_pos = region
    elif isinstance(mut_row, pd.Series):
        chrom = mut_row['Chr']
        pos = mut_row['Pos']
        mut_pos = f"{chrom}:{pos}-{pos} "
        print(mut_pos)
    else:
        mut_pos = ''
    cmd = f"{tool} view {bam_file} {mut_pos} | bam2csv | {HDR_config['editbamdf']} {HDR_config['MINq']}"
    # show_command(cmd)
    bam_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return bam_df


def get_base(read, mut_row, min_q=25):
    '''
    get bases at row position
    '''
    if read['Chr'] != mut_row['Chr']:
        return None
    Seq_pos = mut_row['Start'] - read['Pos'] + read['Soft_start']
    base = read['Seq'][Seq_pos]
    qual = ord(read['Qual'][Seq_pos]) - 33

    if qual >= min_q:
        return 1 if (base == mut_row['Alt']) else 0
    return -1


def get_adjacent_HDR(mut_row, hotspot_df, padding=150):
    '''
    get the adjacent HDR-lanes for each mutation as a HDR_df dataframe for further computation
    '''

    chrom = mut_row['Chr']
    mut_pos = mut_row['Start']
    HDR_df = hotspot_df.query(
        '(Chr == @chrom) and (@mut_pos - @padding < Pos) and (Pos < @mut_pos + @padding)')
    # get distance to mut_spot
    HDR_df.loc[:, 'distance'] = HDR_df['Pos'] - mut_pos
    HDR_df.loc[:, 'Alt'] = HDR_df[['A', 'C', 'T', 'G']].idxmax(axis=1)
    return HDR_df.iloc[:, [0, 1, 2, -1, 3, 10, -2]]


def get_HDR_count(row, df, padding=100):
    HDR_lanes = get_adjacent_HDR(row, df, padding=padding)
    return len(HDR_lanes.index)


def get_intersect_bam(HDR_row, mut_bam, min_q=25):
    '''
    get the reads covering both the mutation and the specific HDR_lane --> intersect_bam
    '''
    pos = HDR_row['Start']
    intersect_bam = mut_bam.query('Pos < @pos < Pos + read_len - Soft_start')
    # prevent key error in next step
    if intersect_bam.empty:
        return intersect_bam
    intersect_bam.loc[:, 'HDRAlt'] = intersect_bam.apply(
        get_base, axis=1, args=(HDR_row,), min_q=min_q)
    # reduce intersect_bam to the reads above threshold quality at that position
    intersect_bam = intersect_bam.query('mutAlt != -1 and HDRAlt != -1')
    return intersect_bam


def get_HDR_base(intersect_read, HDR_row):
    '''
    get bases at mut_row position
    '''
    if intersect_read['Chr'] != HDR_row['Chr']:
        return None
    Seq_pos = HDR_row['Pos'] - intersect_read['Pos'] + \
        intersect_read['Soft_start']
    base = intersect_read['Seq'][Seq_pos]
    qual = ord(intersect_read['Qual'][Seq_pos]) - 33
    if qual > 24:
        # return 1 if this position has an Alt-base
        return int(base == HDR_row['Alt'])
    return 0


def get_covering_reads(bam_df, mut_row):
    '''

    '''
    # get values from mut_row
    pos = mut_row['Start']

    # get reads covering the mutation
    cover_bam = bam_df.query('Pos < @pos < Pos + read_len - Soft_start')
    if len(cover_bam.index) == 0:
        return pd.DataFrame()
    # get reads that are mutated at position
    cover_bam.loc[:, 'mutAlt'] = cover_bam.apply(
        get_base, axis=1, args=(mut_row,))
    return cover_bam


def compute_similarity(HDR_row, cover_bam, min_q=25):
    '''
    for each HDR_row, get the intersecting bam
    '''
    # get the reads covering both the mutation position and the specific HDR_lane --> intersect_bam
    intersect_bam = get_intersect_bam(HDR_row, cover_bam, min_q=min_q)
    if intersect_bam.empty:
        return pd.Series([0, 0, 0, 0, 0], index=['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support'])
    ref_sim = intersect_bam.query('mutAlt == 0 and HDRAlt == 0')[
        'mutAlt'].count()
    ref_support = intersect_bam.query('mutAlt == 0')['mutAlt'].count()
    alt_sim = intersect_bam.query('mutAlt == 1 and HDRAlt == 1')[
        'mutAlt'].count()
    alt_support = intersect_bam.query('mutAlt == 1')['mutAlt'].count()

    # if there is no alt_support, make it zero --> alt_sim/alt_support will be zero as well
    ref_support = ref_support if ref_support else 1
    alt_support = alt_support if alt_support else 1
    support = len(intersect_bam.index)
    result = pd.Series([round(ref_sim / ref_support, 2), ref_support, round(alt_sim / alt_support, 2),
                        alt_support, support], index=['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support'])
    return result


def concat(row):
    ref_support = int(row['RefSupport'])
    ref_sim = int(row['RefSim'] * 100)
    alt_support = int(row['AltSupport'])
    alt_sim = int(row['AltSim'] * 100)
    return f"∆{row['distance']}<Ref:{ref_sim}%({ref_support})><Alt:{alt_sim}%({alt_support})>"


def condense_HDR_info(HDR_df, MinSim=0.9, MinAltSupport=13):
    '''
    reduces the entire HDR_df to entries:
    HDRcount: the number of relevant (similar) lanes around mutation
    HDRmeanSimilarity: the average similarity of these lanes
    HDRinfo: concated string info of all relevant lanes
    '''
    # show_output(HDR_df)
    # select the relevant HDR-lanes / exclude the mutation itself

    HDR_select = HDR_df.query(
        '(AltSupport > @MinAltSupport) and (RefSupport == 0 or RefSim >= @MinSim) and AltSim >= @MinSim')
    if HDR_select.empty:
        return pd.Series([0, 'no similarity in HDR-pattern'], index=['HDRcount', 'HDRinfo'])
    # add info field
    HDR_select.loc[:, 'info'] = HDR_select.apply(concat, axis=1)
    count = HDR_select['info'].count()
    info = HDR_select['info'].str.cat(sep=' | ')
    return pd.Series([count, info], index=['HDRcount', 'HDRinfo'])


def get_HDR_info(mut_row, hotspot_df, bam_df, HDR_config):
    '''
    compute the HDR_info for each mut_row
    --> to be used in filter_HDR.apply
    '''

    show_output(
        f"Analysing Mutation {mut_row['Chr']}:{mut_row['Start']}", multi=True)

    # reduce bam_df to reads covering mutation
    cover_bam = get_covering_reads(bam_df, mut_row)
    # in case there is no coverage on the bam file (should only happen if mutation file is less stringently filtered than filterbam)
    if cover_bam.empty:
        s = pd.Series([0, 'no bam coverage for that mutation'],
                      index=['HDRcount', 'HDRinfo'])
        return s
    # get the HDR_df of adjacent HDR-lanes
    HDR_df = get_adjacent_HDR(mut_row, hotspot_df)
    HDR_df = HDR_df.rename(columns={'Pos': 'Start'})
    # compute the similarity for each HDR_lane
    HDR_df = HDR_df.query('distance != 0')
    if HDR_df.empty:
        s = pd.Series([0, 'no HDR in vincinity'],
                      index=['HDRcount', 'HDRinfo'])
        return s
    HDR_df[['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support']] = HDR_df.apply(
        compute_similarity, axis=1, args=(cover_bam,), min_q=HDR_config['MINq']).fillna(0)
    # HDR_series with fields ['HDRcount', 'HDRmeanSimilarity', 'HDRinfo']
    HDR_series = condense_HDR_info(
        HDR_df, MinSim=HDR_config['MINSIM'], MinAltSupport=HDR_config['MinAltSupport'])
    return HDR_series


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
    filter_pool = Pool(threads)
    # split
    # minimal length of 200 lines
    split_factor = min(math.ceil(len(HDR_df.index) / 200), threads)
    show_output(f"Using {split_factor} cores for filtering HDR-rich mutations")
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
        f'Loaded the bam_file {bam_file} on region {region} for read analysis - {len(bam_df.index)} reads', multi=True)

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


def get_hotspot_df(pileup_file):
    '''
    here comes the pileup computation straight from the bam file
    '''
    pileup_df = get_count_pileup(pileup_file)
    show_output(f"Loading pileup file {pileup_file} finished.")
    hotspot_df = filter_hotspots(pileup_df)
    return hotspot_df
