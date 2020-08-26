import os
import re
import pandas as pd


def sort_df(df):
    '''
    helper for sorting dfs for chromosomes
    '''
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ['chrX', 'chrY']
    df['Chr'] = pd.Categorical(df['Chr'], chrom_list)
    return df.sort_values(['Chr', 'Start'])


def main(s):
    input = str(s.input)
    output = str(s.output)
    params = s.params
    config = s.config

    # ############## COSMIC70 ##############
    cosmic70_dict = {
        'haematopoietic_and_lymphoid_tissue': 5,
        'bone': 2
    }

    # ############## COSMIC91 ##############
    cosmic91_type_score = {
        'acute_myeloid_leukaemia': 8,
        'lymphoid_neoplasm': 2,
        'diffuse_large_B_cell_lymphoma': 3,
        'acute_lymphoblastic_leukaemia': 6,
        'Burkitt_lymphoma': 3,
        'NK-T_cell_lymphoma': 3,
        'chronic_lymphocytic_leukaemia-small_lymphocytic_lymphoma': 3,
        'chronic_myelomonocytic_leukaemia': 6,
        'Hodgkin_lymphoma': 2,
        'acute_leukaemia_of_ambiguous_lineage': 5,
        'acute_lymphoblastic_B_cell_leukaemia': 3,
        'acute_lymphoblastic_T_cell_leukaemia': 3,
        'acute_myeloid_leukaemia_associated_with_MDS': 6,
        'blastic_plasmacytoid_dendritic_cell_neoplasm': 3,
        'chronic_myeloid_leukaemia': 6,
        'juvenile_myelomonocytic_leukaemia': 6,
        'myelodysplastic_syndrome': 5,
        'plasma_cell_myeloma': 3,
        'acute_myeloid_leukaemia_therapy_related': 6,
        'adult_T_cell_lymphoma-leukaemia': 4,
        'blast_phase_chronic_myeloid_leukaemia': 6,
        'mast_cell_neoplasm': 6,
        'follicular_lymphoma': 3,
        'T_cell_large_granular_lymphocytic_leukaemia': 3,
        'marginal_zone_lymphoma': 3,
        'benign': -1,
        'normal': -1
    }
    cosmic91_location_score = {
        'haematopoietic_and_lymphoid_tissue': 3,
        'femur': 3,
        'bone': 3,
        'thymus': 2,
        'tonsil': 2,
        'tibia': 2,
        'medulla': 3,  # is this the bone marrow
        'spleen': 2
    }

    # ############## CLINVAR #################
    CLNDN_factorial = {
        'carcinoma': 2,
        'neoplasm': 2,
        'cancer': 2,
        'malignant': 1.5,
        'melanoma': 1.5,
        'myeloma': 5,
        'lymphoma': 5,
        'lymphatic': 3,
        'blastoma': 1.5,
        'immunodeficiency': 1.25
    }

    CLNSIG_score = {
        'Affects': 0.2,
        'Benign': -1,
        'Benign/Likely_benign': -0.75,
        'Conflicting_interpretations_of_pathogenicity': 0,
        'Likely_benign': -0.5,
        'Likely_pathogenic': 0.5,
        'Pathogenic': 1,
        'Pathogenic/Likely_pathogenic': 0.75,
        'Uncertain_significance': 0
    }

    # ############## ====> CLINSCORE #################
    ClinScore = {
        'cosmic70_score': 2,  # derived score
        'cosmic91_score': 2,  # derived score
        'clinvar_score': 2,  # derived score
        'icgc29_freq': 2000  # derived score --> adjust
    }

    # import the scripts
    extended_output = params.extended_output
    candidate_list = params.candidate_list
    driver_list = params.driver_list
    hotspot_list = params.hotspot_list

    print(f'Started editing and basic filtering for {input}.')
    anno_df = pd.read_csv(input, sep='\t')

    ############# PREDICTIONS TO KEEP ################
    predictions = ["Polyphen2", "SIFT", "MutationTaster"]

    #####################################################################
    # #################### CUSTOM WEIGHTS ################################
    '''
    impacts how the cosmic string is translated into scalar value
    '''

    # only get the non-zero keys
    clinscore_cols = [col for col in ClinScore.keys() if ClinScore[col]]

    #####################################################################

    # ########## CLINICAL ROWS ##############
    # def is_clin_col(col):
    #     for key in ['cosmic', 'CLN', 'Clin', 'icgc']:
    #         if key in col:
    #             return True
    #     return False

    def get_PoN_info(df):
        print('Calculating PoN metrix')
        df.loc[:, 'PoN-Alt-Sum'] = df['PoN-Alt'].str.replace('-', '|').str.split(
            '|').apply(lambda array: sum([int(count) for count in array]))
        df.loc[:, 'PoN-Ref-Sum'] = df['PoN-Ref'].str.replace('-', '|').str.split(
            '|').apply(lambda array: sum([int(count) for count in array]))
        df.loc[:, 'PoN-Alt-NonZeros'] = df['PoN-Alt'].str.count(
            r'[0-9]+') - df['PoN-Alt'].str.count('0')
        df.loc[:, 'PoN-Ratio'] = df['PoN-Alt-Sum'] / df['PoN-Ref-Sum']
        return df

    def resort_cols(df):
        '''
        resort the columns and removes prediction columns based on
        '''
        cols = list(df.columns)

        # #### DEBUGGING ######
        # print('All Columns:')
        # for i, col in enumerate(cols):
        #     print(i, col)
        ######################
        start_cols = cols[:11]
        quant_cols = cols[11:26] + ['FisherScore', 'EBscore', 'PoN-Ref',
                                    'PoN-Alt', 'PoN-Ref-Sum', 'PoN-Alt-Sum', 'PoN-Alt-NonZeros', 'PoN-Ratio']
        if 'A|a|G|g|C|c|T|t|I|i|D|d' in cols:
            quant_cols.append('A|a|G|g|C|c|T|t|I|i|D|d')
        clin_cols = ['ClinScore', 'cosmic91_ID', 'cosmic91_type', 'cosmic91_score', 'cosmic70_ID', 'cosmic70_freq',
                     'cosmic70_type', 'cosmic70_score', 'CLNALLELEID', 'CLNDN', 'CLNSIG', 'clinvar_score', 'icgc29_ID', 'icgc29_freq']
        pop_cols = [col for col in cols if re.match(
            "|".join(["avsnp", "dbSNP", "esp6500", "1000g", "gnomAD"]), col)]
        print(f"Using population data: {' '.join(pop_cols)}")
        # keep only the pred_cols defined in predictions list above
        pred_cols = [col for col in cols if re.match(
            "|".join(predictions), col)]

        pred_cols = cols[47:-13] if extended_output else pred_cols
        # 13 <== the added extracted and score columns make up 8 columns + 4 PoN columns = 12:
        # 4:    'icgc29_freq'
        # 5-7:  'clinvar_score', 'cosmic70_score', 'cosmic91_score'
        # 8:    'ClinScore'
        # 9-12: PoN-info (4 columns)
        print(f"Keeping predictions: {' '.join(pred_cols)}")
        new_cols = start_cols + quant_cols + clin_cols + pop_cols + pred_cols
        return df[new_cols]

    def get_clinical_scores(df):
        '''
        extract, score and realign clinical columns
        '''
        # ############## CONVERSION OF COLUMNS ##############################
        # ############## ICGC29 ###########
        def addICGC(df):
            ICGC = df['icgc29_Affected'].str.extract(
                r'^([0-9]+)/([0-9]+)$').fillna(0).astype('int')
            df['icgc29_freq'] = (ICGC[0] / ICGC[1]).round(4).fillna(".")
            return df.drop(columns='icgc29_Affected')

        # ############## COSMIC70 #########
        def addCosmic70(df):
            pattern = r'(?:ID=(?P<cosmID>COSM[0-9]+(?:,COSM[0-9]+)?);OCCURENCE=)?(?P<freq>[0-9]+)\((?P<organ>[A-Z_a-z]+)\)'
            df[['cosmic70_ID', 'cosmic70_freq', 'cosmic70_type']] = df['cosmic70'].str.extractall(pattern).astype({'cosmID': 'str', 'freq': 'int', 'organ': 'str'}).reset_index(
                'match').drop(columns='match').reset_index().groupby('index').aggregate({'cosmID': 'min', 'freq': 'sum', 'organ': lambda col: col.str.cat(sep='+')})
            df.loc[:, 'cosmic70_freq'] = df['cosmic70_freq'].fillna(
                0).astype('int')
            df.loc[:, 'cosmic70_ID'] = df['cosmic70_ID'].fillna('.')
            df.loc[:, 'cosmic70_type'] = df['cosmic70_type'].fillna('.')
            return df.drop(columns='cosmic70')

        # ############## COMPUTATION OF SCORES ##############################
        def cosmic70score(row):
            if row['cosmic70_type'] != ".":
                score = 1
                for location in cosmic70_dict.keys():
                    if location in row['cosmic70_type']:
                        score += cosmic70_dict[location]
                return row['cosmic70_freq'] * score
            else:
                return 0

        def cosmic91_score(row):
            return (1 + cosmic91_type_score.get(row['types'], 0) + cosmic91_location_score.get(row['location'], 0)) * (row['types'] != ".") * int(row['count'])

        cosmic91_pattern = r'(?P<count>[0-9]+)x\((?P<types>[^0-9@)]+)@(?P<location>[^0-9@)]+)\)'

        def get_CLINVARscore(row):
            '''
            converts the CLINVAR info into scalar clinvar_score
            '''

            def CLNDN2score(clndn):
                '''
                accumulates a factor for multiplication with CLNSIG_score
                '''

                factor = 1
                for key in CLNDN_factorial:
                    if key in clndn:
                        factor *= CLNDN_factorial[key]
                return factor

            if row['CLNDN'] == ".":
                return 0
            return 1 + CLNDN2score(row['CLNDN']) + CLNSIG_score.get(row['CLNSIG'].split(',')[0], 0)

        # INFERRED COLUMNS
        print('Add Cosmic70 derived columns')
        df = addCosmic70(df)
        print('Add ICGC derived columns')
        df = addICGC(df)

        print('Derive Clinical Scores from Cosmic and Clinvar')
        # SCALAR SCORES FROM CLINICAL DBS
        df['clinvar_score'] = df.apply(get_CLINVARscore, axis=1)
        df['cosmic70_score'] = df.apply(cosmic70score, axis=1)
        df['cosmic91_score'] = df['cosmic91_type'].str.extractall(cosmic91_pattern).apply(
            cosmic91_score, axis=1).reset_index().drop(columns='match').groupby('level_0').sum()
        df['cosmic91_score'] = df['cosmic91_score'].fillna(0)

        # GET COMBINED ClinScore
        df['ClinScore'] = 0
        print('      Combining clinical scores into ClinScore')
        for col in clinscore_cols:
            print("            ", col)
            df.loc[df[col] != ".", 'ClinScore'] += ClinScore[col] * \
                df[df[col] != "."][col]
        return df

    def get_gene_lists(df, candidate_list, driver_list, hotspot_list):
        '''
        add gene info from attached lists
        '''

        # ## candidate list #########
        list_cols = []
        STATIC = config['paths']['mystatic']

        if candidate_list:
            candidate_file = os.path.join(STATIC, candidate_list)
            print(f'Using candidate list {candidate_file}')
            candidate_list = list(pd.read_csv(candidate_file, header=None)[0])
            df['isCandidate'] = df['Gene'].isin(candidate_list).astype(int)
            list_cols.append('isCandidate')
        if driver_list:
            driver_file = os.path.join(STATIC, driver_list)
            print(f'Using candidate list {driver_file}')
            driver_list = list(pd.read_csv(driver_file, header=None)[0])
            df['isDriver'] = df['Gene'].isin(driver_list).astype(int)
            list_cols.append('isDriver')

        # add the hotspot mutations to the file
        if hotspot_list:
            hotspot_file = os.path.join(STATIC, hotspot_list)
            hotspots = pd.read_csv(hotspot_file, sep='\t').drop(
                columns=['Protein'])
            df = df.merge(hotspots, how='left', on=[
                          'Chr', 'Start', 'Ref', 'Alt', 'Gene'])
            list_cols += ['ChipID', 'ChipPub', 'ChipFreq']

        # resort columns
        cols = list(df.columns)
        start_cols = cols[:11]
        # the last columns are the newly addded columns isDriver and isCandidate
        # have to be omitted
        rest_cols = cols[11:-len(list_cols)]
        new_cols = start_cols + list_cols + rest_cols
        return df[new_cols]

    pon_info_df = get_PoN_info(anno_df)
    clin_df = get_clinical_scores(pon_info_df)
    # RESORT THE COLUMNS
    sorted_df = resort_cols(clin_df)
    candidate_df = get_gene_lists(
        sorted_df, candidate_list, driver_list, hotspot_list)

    # sort by chrom
    candidate_df = sort_df(candidate_df)
    # this is raw unfiltered data, only informative columns were added
    candidate_df.to_csv(output, sep='\t', index=False)
    print(f"Writing edited mutation list with added columns to {output}.")


if __name__ == "__main__":
    main(snakemake)
