import pandas as pd
import pysam
from script_utils import show_output


def get_count_pileup(filename, _type='Tumor'):
    cols = [0, 1, 2, 3, 4, 5] if _type == 'Tumor' else [0, 1, 2, 6, 7, 8]
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
    hotspot_df = pileup_df.query('(AltSum >= @minAlt) and (@minRatio <= AltRatio <= @maxRatio)')
    return hotspot_df

# cigar_pattern = re.compile(r'^([0-9]+)([NMDIS])(.*)')


def bam2df(bam_file, q=20):
    '''
    reads the sub bam into a df using pysam
    pysam is only used to extract the read columns in a sensible way
    '''

    lst = []
    # read the bam file line by line into a list of dicts
    with pysam.AlignmentFile(bam_file, "r") as bam_file:
        for i, line in enumerate(bam_file):
            row = line.to_dict()
            # extract the tag elements into dictionary keys
            row.update({tag.split(':')[0]: tag.split(':')[2]
                        for tag in row['tags']})
            row.pop('tags')
            lst.append(row)
    # read list into a dataframe
    bam_df = pd.DataFrame(lst)
    bam_df['map_quality'] = bam_df['map_quality'].astype(int)
    return bam_df.query('map_quality >= @q')


def editbamdf(df):
    '''
    clean_up the bam file
    '''

    # convert position into integer
    df['ref_pos'] = df['ref_pos'].astype(int)
    df['read_len'] = df['seq'].str.len()  # only for debugging
    # extract the intron sizes from the cigar string
    df['Chr_len'] = df['read_len']
    # extraction of soft-clipped bases
    df['soft_start'] = df['cigar'].str.extract(
        r'(^[0-9]+)S').fillna(0).astype(int)
    # filter out reads without a Cell barcode and return only desired columns
    df = df[['name', 'ref_name', 'ref_pos', 'read_len',
             'seq', 'qual', 'soft_start', 'cigar']]
    df.columns = ['name', 'Chr', 'Pos', 'read_len',
                  'Seq', 'Qual', 'Soft_start', 'Cigar']
    return df


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
    HDR_df = hotspot_df.query('(Chr == @chrom) and (@mut_pos - @padding < Pos) and (Pos < @mut_pos + @padding)')
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


def condense_HDR_info(HDR_df, MINSIM=0.9):
    '''
    reduces the entire HDR_df to entries:
    HDRcount: the number of relevant (similar) lanes around mutation
    HDRmeanSimilarity: the average similarity of these lanes
    HDRinfo: concated string info of all relevant lanes
    '''
    # show_output(HDR_df)
    # select the relevant HDR-lanes / exclude the mutation itself

    HDR_select = HDR_df.query(
        'AltSupport > 13 and (RefSupport == 0 or RefSim >= @MINSIM) and AltSim >= @MINSIM')
    if HDR_select.empty:
        return pd.Series([0, 'no similarity in HDR-pattern'], index=['HDRcount', 'HDRinfo'])
    # add info field
    HDR_select.loc[:, 'info'] = HDR_select.apply(concat, axis=1)
    count = HDR_select['info'].count()
    info = HDR_select['info'].str.cat(sep=' | ')
    return pd.Series([count, info], index=['HDRcount', 'HDRinfo'])


def get_HDR_info(mut_row, hotspot_df, bam_df, MINSIM=0.9, min_q=25):
    '''
    compute the HDR_info for each mut_row
    --> to be used in filter_HDR.apply
    '''
    show_output(f"Analysing Mutation {mut_row['Chr']}:{mut_row['Start']}")

    # reduce bam_df to reads covering mutation
    cover_bam = get_covering_reads(bam_df, mut_row)
    # in case there is no coverage on the bam file (should only happen if mutation file is less stringenty filtered than filterbam)
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
        compute_similarity, axis=1, args=(cover_bam,), min_q=min_q).fillna(0)
    # HDR_series with fields ['HDRcount', 'HDRmeanSimilarity', 'HDRinfo']
    HDR_series = condense_HDR_info(HDR_df, MINSIM=MINSIM)
    return HDR_series


def get_HDR(bam_file, mut_df, pileup_file, _type, MINSIM, padding, min_HDR_count, MINQ):

    # get the right pileup_df (different cols for Tumor and Normal)
    pileup_df = get_count_pileup(pileup_file, _type=_type)
    show_output(f"Loading pileup file {pileup_file} finished.")
    hotspot_df = filter_hotspots(pileup_df)
    show_output(
        f"Detected {len(hotspot_df.index)} putative HDR lanes in {bam_file}.")
    # enumerate the HDRs in vicinity of mutations
    mut_df.loc[:, 'HDR'] = mut_df.apply(
        get_HDR_count, axis=1, args=(hotspot_df,), padding=padding)

    # continue with HDR-rich mutations
    filter_HDR = mut_df.query('HDR >= @min_HDR_count')
    show_output(f"Found {len(filter_HDR.index)} HDR-rich mutations")

    # ###### BAM ANALYSIS ##################################
    # get the bam_df for analysis in single-read resolution
    bam_df = editbamdf(bam2df(bam_file))
    show_output(f'Loaded the bam_file {bam_file} for read analysis')
    mut_df[['HDRcount', 'HDRinfo']] = filter_HDR.apply(
        get_HDR_info, axis=1, args=(hotspot_df, bam_df), MINSIM=MINSIM, min_q=MINQ)
    mut_df.loc[:, 'HDRcount'] = mut_df['HDRcount'].fillna(0).astype(int)
    mut_df.loc[:, 'HDRinfo'] = mut_df['HDRinfo'].fillna('no HDR')
    return mut_df
