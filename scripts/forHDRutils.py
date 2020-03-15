

def get_covering_reads(bam_df, mut_row):
    '''

    '''
    # get values from mut_row
    pos = mut_row['Start']
    Alt = mut_row['Alt']

    # get reads covering the mutation
    cover_bam = bam_df.query('Pos < @pos < Pos + read_len - Soft_start')
    if len(cover_bam.index) == 0:
        return pd.DataFrame()
    # get reads that are mutated at position
    cover_bam['mutAlt'] = cover_bam.apply(get_base, axis=1, args=(mut_row,))
    return cover_bam


def compute_similarity(HDR_row, cover_bam):
    '''
    for each HDR_row, get the intersecting bam
    '''
    # get the reads covering both the mutation position and the specific HDR_lane --> intersect_bam
    intersect_bam = get_intersect_bam(HDR_row, cover_bam)
    if intersect_bam.empty:
        return pd.Series([0, 0, 0, 0, 0], index=['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support'])
    ref_sim = intersect_bam.query('mutAlt == 0 and HDRAlt == 0')[
        'mutAlt'].count()
    ref_support = intersect_bam.query('mutAlt == 0')['mutAlt'].count()
    alt_sim = intersect_bam.query('mutAlt == 1 and HDRAlt == 1')[
        'mutAlt'].count()
    alt_support = intersect_bam.query('mutAlt == 1')['mutAlt'].count()
    support = len(intersect_bam.index)
    result = pd.Series([round(ref_sim/ref_support, 2), ref_support, round(alt_sim/alt_support, 2),
                        alt_support, support], index=['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support'])
    return result


def concat(row):
    ref_support = int(row['RefSupport'])
    ref_sim = int(row['RefSim'] * 100)
    alt_support = int(row['AltSupport'])
    alt_sim = int(row['AltSim'] * 100)
    return f"∆{row['distance']}<Ref:{ref_sim}%({ref_support})><Alt:{alt_sim}%({alt_support})>"


def condense_HDR_info(HDR_df):
    '''
    reduces the entire HDR_df to entries:
    HDRcount: the number of relevant (similar) lanes around mutation
    HDRmeanSimilarity: the average similarity of these lanes
    HDRinfo: concated string info of all relevant lanes
    '''
    # print(HDR_df)
    # select the relevant HDR-lanes / exclude the mutation itself

    HDR_select = HDR_df.query(
        'AltSupport > 13 and (RefSupport == 0 or RefSim >= @MINSIM) and AltSim >= @MINSIM')
    if HDR_select.empty:
        return pd.Series([0, 'no similarity in HDR-pattern'], index=[f'HDRcount', f'HDRinfo'])
    # add info field
    HDR_select.loc[:, 'info'] = HDR_select.apply(concat, axis=1)
    count = HDR_select['info'].count()
    info = HDR_select['info'].str.cat(sep=' | ')
    return pd.Series([count, info], index=['HDRcount', 'HDRinfo'])


####### get_HDR_info for debugging #########################
def get_HDR_info(mut_row, hotspot_df, bam_df):
    '''
    compute the HDR_info for each mut_row 
    --> to be used in filter_HDR.apply
    '''
    print(f"Analysing Mutation {mut_row['Chr']}:{mut_row['Start']}")

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
                      index=[f'HDRcount', f'HDRinfo'])
        return s
    HDR_df[['RefSim', 'RefSupport', 'AltSim', 'AltSupport', 'support']] = HDR_df.apply(
        compute_similarity, axis=1, args=(cover_bam,)).fillna(0)
    # HDR_series with fields ['HDRcount', 'HDRmeanSimilarity', 'HDRinfo']
    HDR_series = condense_HDR_info(HDR_df)
    return HDR_series


def getHDR(bam_file, filter_df):
    print(f"Piling up {bam_file}..")
    pileup_file = get_clean_pileup(bam_file)
    pileup_df = get_count_pileup(pileup_file)
    print(f"Pileup finished.")
    hotspot_df = filter_hotspots(pileup_df)
    print(
        f"Detected {len(hotspot_df.index)} putative HDR lanes in {bam_file}.")
    # enumerate the HDRs in vicinity of mutations
    filter_df['HDR'] = filter_df.apply(
        get_HDR_count, axis=1, args=(hotspot_df,))

    # continue with HDR-rich mutations
    filter_HDR = filter_df.query('HDR > 0')
    print(f"Found {len(filter_HDR.index)} HDR-rich mutations")

    ####### BAM ANALYSIS ##################################
    # get the bam_df for analysis in single-read resolution
    bam_df = editbamdf(bam2df(bam_file))
    print(f'Loaded the bam_file {bam_file} for read analysis')
    filter_df[['HDRcount', 'HDRinfo']] = filter_HDR.apply(
        get_HDR_info, axis=1, args=(hotspot_df, bam_df))
    filter_df['HDRcount'] = filter_df['HDRcount'].fillna(0).astype(int)
    filter_df['HDRinfo'] = filter_df['HDRinfo'].fillna('no HDR')
    return filter_df


def masterHDR(mut_df, tumor_bam='', normal_bam=''):

    bam_dict = {'Tumor': tumor_bam, 'Normal': normal_bam}
    ####### PILEUP ANALYSIS ##############################
    for type in ['Tumor', 'Normal']:
        print(f'Analysing {type}')
        mut_df = getHDR(bam_dict[type], mut_df)
        mut_df = mut_df.rename(columns={
                               'HDR': f'{type}HDRcand', 'HDRcount': f'{type}HDRcount', 'HDRinfo': f'{type}HDRinfo'})
    return mut_df
