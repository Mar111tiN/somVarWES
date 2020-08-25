import pandas as pd
from script_utils import show_output


def get_base(read, mut_row, min_q=25):
    '''
    get bases at row position
    '''
    # chrom check is not neccessary
    # if read['Chr'] != mut_row['Chr']:
    #     return None
    Seq_pos = mut_row['Start'] - read['Pos'] + read['Soft_start']
    base = read['Seq'][Seq_pos]
    qual = ord(read['Qual'][Seq_pos]) - 33

    if qual >= min_q:
        return 1 if (base == mut_row['Alt']) else 0
    return -1


def get_covering_reads(bam_df, mut_row, HDR_config):
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
        get_base, axis=1, args=(mut_row, ), min_q=HDR_config['MINQ'])
    # !!!!!!!!!!!
    # should I reduce already here to good reads
    cover_bam = cover_bam.query('mutAlt != -1')
    return cover_bam


def get_adjacent_HDR2(mut_row, hotspot_df, padding=150):
    '''
    get the adjacent HDR-lanes for each mutation as a HDR_df dataframe for further computation
    '''

    chrom = mut_row['Chr']
    mut_pos = mut_row['Start']
    HDR_df = hotspot_df.query(
        '(@mut_pos - @padding < Pos < @mut_pos + @padding)')
    # get Dist to mut_spot
    HDR_df.loc[:, 'Dist'] = HDR_df['Pos'] - mut_pos
    cols = ['Chr', 'Pos', 'Ref', 'Depth', 'Alt', 'AltSum', 'Dist']
    return HDR_df[cols]


def get_HDR_count(row, df, padding=100):
    HDR_lanes = get_adjacent_HDR2(row, df, padding=padding)
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
    return f"∆{row['Dist']}<Ref:{ref_sim}%({ref_support})><Alt:{alt_sim}%({alt_support})>"


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
    cover_bam = get_covering_reads(bam_df, mut_row, HDR_config=HDR_config)
    # in case there is no coverage on the bam file (should only happen if mutation file is less stringently filtered than filterbam)
    if cover_bam.empty:
        s = pd.Series([0, 'no bam coverage for that mutation'],
                      index=['HDRcount', 'HDRinfo'])
        return s
    # get the HDR_df of adjacent HDR-lanes
    HDR_df = get_adjacent_HDR2(
        mut_row, hotspot_df, padding=HDR_config['PAD']+10)
    HDR_df = HDR_df.rename(columns={'Pos': 'Start'})
    # compute the similarity for each HDR_lane
    HDR_df = HDR_df.query('Dist != 0')
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
