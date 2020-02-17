import pandas as pd
import pysam
from script_utils import show_output


def get_count_pileup(filename):
    pileup = pd.read_csv(
        filename,
        sep='\t',
        header=None,
        usecols=[0, 1, 2, 3, 4, 5],
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
            row.update({tag.split(':')[0]: tag.split(':')[2] for tag in row['tags']})
            row.pop('tags')
            lst.append(row)
    # read list into a dataframe
    bam_df = pd.DataFrame(lst)
    bam_df['map_quality'] = bam_df['map_quality'].astype(int)
    return bam_df.query('map_quality >= @q')


def editbamdf(df):
    '''
    reads the sub bam into a df using pysam
    pysam is only used to extract the read columns in a sensible way
    '''

    # convert position into integer
    df['ref_pos'] = df['ref_pos'].astype(int)
    df['read_len'] = df['seq'].str.len()  # only for debugging
    # extract the intron sizes from the cigar string
    df['Chr_len'] = df['read_len']
    # extraction of soft-clipped bases
    df['soft_start'] = df['cigar'].str.extract(r'(^[0-9]+)S').fillna(0).astype(int)
    # filter out reads without a Cell barcode and return only desired columns
    df = df[['name', 'ref_name', 'ref_pos', 'read_len', 'seq', 'qual', 'soft_start', 'cigar']]
    df.columns = ['name', 'Chr', 'Pos', 'read_len', 'Seq', 'Qual', 'Soft_start', 'Cigar']
    return df


def get_base(read, mut_row):
    '''
    get bases at mut_row position
    '''
    if read['Chr'] != mut_row['Chr']:
        return None
    Seq_pos = mut_row['Start'] - read['Pos'] + read['Soft_start']
    base = read['Seq'][Seq_pos] 
    qual = ord(read['Qual'][Seq_pos]) - 33
    if qual > 24:
        return base
    return "N"


def get_mutated_bam(bam_df, mut_row):
    # get values from mut_row
    pos = mut_row['Start']
    Alt = mut_row['Alt']

    # get reads covering the mutation
    cover_bam = bam_df.query('Pos < @pos < Pos + read_len - Soft_start')

    # get reads that are mutated at position
    cover_bam.loc[:, 'mut_pos'] = cover_bam.apply(get_base, axis=1, args=(mut_row,))
    mut_bam = cover_bam.query('mut_pos == @Alt')
    return mut_bam


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


def get_intersect_bam(HDR_row, mut_bam):
    '''
    get the reads covering both the mutation and the specific HDR_lane --> intersect_bam
    '''
    pos = HDR_row['Pos']
    chrom = HDR_row['Chr']
    Alt = HDR_row['Alt']
    intersect_bam = mut_bam.query('Pos < @pos < Pos + read_len - Soft_start')
    return intersect_bam


def get_HDR_base(intersect_read, HDR_row):
    '''
    get bases at mut_row position
    '''
    if intersect_read['Chr'] != HDR_row['Chr']:
        return None
    Seq_pos = HDR_row['Pos'] - intersect_read['Pos'] + intersect_read['Soft_start']
    base = intersect_read['Seq'][Seq_pos]
    qual = ord(intersect_read['Qual'][Seq_pos]) - 33
    if qual > 24:
        return int(base == HDR_row['Alt']) # return 1 if this position has an Alt-base
    return 0


def get_HDR_Alt_percentage(HDR_row, intersect_bam):
    Alt_sum = intersect_bam.apply(get_HDR_base, axis=1, args=(HDR_row,)).sum()
    HDR_Alt_percentage = Alt_sum / len(intersect_bam.index)
    return HDR_Alt_percentage.round(2)


def compute_similarity(HDR_row, mut_bam):
    '''
    for each HDR_row, get the intersecting bam
    '''
    intersect_bam = get_intersect_bam(HDR_row, mut_bam)
    if intersect_bam.empty:
        return pd.Series([0,0], index=['similarity', 'support'])
    return pd.Series([get_HDR_Alt_percentage(HDR_row, intersect_bam), len(intersect_bam.index)], index=['similarity', 'support'])


def concat(row):
    support = int(row['support'])
    simil = int(row['similarity'] * 100)
    return f"∆{row['distance']}:{simil}%({support})"


def condense_HDR_info(HDR_df, min_sim=0.85):
    '''
    reduces the entire HDR_df to entries:
    HDRcount: the number of relevant (similar) lanes around mutation
    HDRmeanSimilarity: the average similarity of these lanes
    HDRinfo: concated string info of all relevant lanes
    '''
    # print(HDR_df)
    # select the relevant HDR-lanes / exclude the mutation itself
    HDR_select = HDR_df.query('(distance != 0) and (similarity > @min_sim) and (support > 9)')
    if HDR_select.empty:
        return pd.Series([0, 0, 'no similarity in HDR-pattern'], index=[
            'HDRcount',
            'HDRmeanSimilarity',
            'HDRinfo'
            ])
    # add info field
    HDR_select.loc[:, 'info'] = HDR_select.apply(concat, axis=1)
    count = HDR_select['info'].count()
    mean_sim = HDR_select['similarity'].mean()
    info = HDR_select['info'].str.cat(sep=' | ')
    return pd.Series([count, mean_sim, info], index=['HDRcount', 'HDRmeanSimilarity', 'HDRinfo'])


def get_HDR_info(mut_row, hotspot_df, bam_df, padding=100, min_sim=0.85):
    '''
    compute the HDR_info for each mut_row 
    --> to be used in filter_HDR.apply
    '''
    show_output(f"Analysing Mutation {mut_row['Chr']}:{mut_row['Start']}", time=False)
    # reduce bam_df to mutation vicinity
    mut_bam = get_mutated_bam(bam_df, mut_row)
    # get the HDR_df of adjacent HDR-lanes
    HDR_df = get_adjacent_HDR(mut_row, hotspot_df, padding=padding)

    # compute the similarity for each HDR_lane
    HDR_df[['similarity', 'support']] = HDR_df.apply(compute_similarity, axis=1, args=(mut_bam,))
    # HDR_series with fields ['HDRcount', 'HDRmeanSimilarity', 'HDRinfo']
    HDR_series = condense_HDR_info(HDR_df, min_sim=min_sim)
    return HDR_series


def masterHDR(pileup_file='', bam_file='', filter_df=None, min_sim=.90, padding=100):
    '''
    compute and output the HDR for the mutation file
    '''

    # get the HDR_hotspots in the samples bam file
    pileup_df = get_count_pileup(pileup_file)
    show_output(f"Imported pileup {pileup_file} .", time=False)
    hotspot_df = filter_hotspots(pileup_df)
    show_output(f"Detected {len(hotspot_df.index)} putative HDR lanes in {bam_file}.")
    # enumerate the HDRs in vicinity of mutations
    filter_df['HDR'] = filter_df.apply(
        get_HDR_count,
        axis=1,
        args=(hotspot_df,),
        padding=padding
        )

    # continue with HDR-rich mutations
    filter_HDR = filter_df.query('HDR > 1')
    print(f"Found {len(filter_HDR.index)} HDR-rich mutations")
    # get the bam_df for analysis in single-read resolution
    bam_df = editbamdf(bam2df(bam_file))
    print(f'Loaded the bam_file {bam_file} for read analysis')
    filter_df[['HDRcount', 'HDRmeanSimilarity', 'HDRinfo']] = filter_HDR.apply(
        get_HDR_info,
        axis=1,
        args=(hotspot_df, bam_df),
        padding=padding,
        min_sim=min_sim
        )
    filter_df['HDRcount'] = filter_df['HDRcount'].fillna(0).astype(int)
    filter_df['HDRmeanSimilarity'] = filter_df['HDRmeanSimilarity'].fillna(0)
    filter_df['HDRinfo'] = filter_df['HDRinfo'].fillna('no HDR')
    return filter_df
    # return filter_HDR, hotspot_df, bam_df
