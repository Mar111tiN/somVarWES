#!/usr/bin/env python

import os
import argparse
import pandas as pd


# set_up the parser for input
parser = argparse.ArgumentParser('filters annovar output with custom criteria')
parser.add_argument('-TM2_limit', type=int, help='test value for argparse')
parser.add_argument('-static_path', type=str, default='', help='path to the static folder containing the gene lists')
parser.add_argument('-candidate_list', type=str, default='', help='path to a file containing candidate genes')
parser.add_argument('-driver_list', type=str, default='', help='path to a file containing known driver genes')
parser.add_argument('-keep_syn', type=bool, default=False, help='True if you want to keep synonymous mutations')
parser.add_argument('-threads', type=int, default=1, help='threads for assuring memory')
parser.add_argument('input', type=str, help='input from annovar')
parser.add_argument('output', type=str, help='output file to filtered/..')


# read arguments
args = parser.parse_args()
keep_syn = args.keep_syn
i, o = args.input, args.output
threads = args.threads
keep_syn = args.keep_syn
output_base = os.path.splitext(o)[0]

# params
TM2_limit = args.TM2_limit

print(f'Started editing and basic filtering for {i}.')
anno_df = pd.read_csv(i, sep='\t')

#####################################################################
# ############## ICGC29 ##############################################


def addICGC(df):
    ICGC = df['icgc29_Affected'].str.extract(r'^([0-9]+)/([0-9]+)$').fillna(0).astype('int')
    df['icgc29_freq'] = (ICGC[0] / ICGC[1]).fillna(".")
    return df.drop(columns='icgc29_Affected')


#####################################################################
# ############## COSMIC70 ############################################

def addCosmic70(df):
    pattern = r'(?:ID=(?P<cosmID>COSM[0-9]+(?:,COSM[0-9]+)?);OCCURENCE=)?(?P<freq>[0-9]+)\((?P<organ>[A-Z_a-z]+)\)'
    df[['cosmic70_ID', 'cosmic70_freq', 'cosmic70_type']] = df['cosmic70'].str.extractall(pattern).astype({'cosmID':'str', 'freq':'int', 'organ':'str'}).reset_index('match').drop(columns='match').reset_index().groupby('index').aggregate({'cosmID':'min', 'freq':'sum', 'organ': lambda col: col.str.cat(sep='+')})
    df.loc[:, 'cosmic70_freq'] = df['cosmic70_freq'].fillna(0).astype('int')
    df.loc[:, 'cosmic70_ID'] = df['cosmic70_ID'].fillna('.')
    df.loc[:, 'cosmic70_type'] = df['cosmic70_type'].fillna('.')
    return df.drop(columns='cosmic70')


cosmic70_dict = {
    'haematopoietic_and_lymphoid_tissue': 5,
    'bone': 2
}


def cosmic70score(row):
    if row['cosmic70_type'] != ".":
        score = 1
        for location in cosmic70_dict.keys():
            if location in row['cosmic70_type']:
                score += cosmic70_dict[location]
        return row['cosmic70_freq'] * score
    else:
        return 0


#####################################################################
############### COSMIC90 ############################################

cosmic90_type_score = {
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

cosmic90_location_score = {
    'haematopoietic_and_lymphoid_tissue': 3,
    'femur': 3,
    'bone': 3,
    'thymus': 2,
    'tonsil': 2,
    'tibia': 2,
    'medulla': 3,  # is this the bone marrow
    'spleen': 2
}

cosmic90_score = lambda row: (1 + cosmic90_type_score.get(row['types'], 0) + cosmic90_location_score.get(row['location'], 0)) * (row['types'] != ".") * int(row['count'])
pattern = r'(?P<count>[0-9]+)x\((?P<types>[^0-9@)]+)@(?P<location>[^0-9@)]+)\)'


#####################################################################
# ############## CLINVAR #############################################

CLNDN_factorial = {
    'carcinoma': 2,
    'neoplasm': 2,
    'cancer': 2,
    'malignant': 1.5,
    'melanoma': 1.5,
    'myeloma': 5,
    'lymphoma': 5,
    'lymphatic':3,
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


def CLNDN2score(clndn):
    '''
    accumulates a factor for multiplication with CLNSIG_score
    '''

    factor = 1
    for key in CLNDN_factorial:
        if key in clndn:
            factor *= CLNDN_factorial[key]
    return factor


def get_CLINVARscore(row):
    if row['CLNDN'] == ".":
        return 0
    return 1 + CLNDN2score(row['CLNDN']) + CLNSIG_score.get(row['CLNSIG'].split(',')[0], 0)


#####################################################################
# ############## COMBINED CLINSCORE ##################################

# ########## CUSTOM WEIGHTS ##############

# derived clinical significance scores --> ClinScore
ClinScore = {
    'cosmic70_score': 2,  # derived score
    'cosmic90_score': 2,  # derived score
    'clinvar_score': 2,  # derived score
    'icgc29_freq': 200  # derived score --> should be much higher!!!
}

# only get the non-zero keys
clinscore_cols = [col for col in ClinScore.keys() if ClinScore[col]]


#####################################################################
# ############## GET ALL CLINICAL SCORES #############################

# ########## CLINICAL ROWS ##############
def is_clin_col(col):
    for key in ['cosmic', 'CLN', 'Clin', 'icgc']:
        if key in col:
            return True
    return False


def resort_cols(df):
    # resort the columns
    cols = list(df.columns)
    start_cols = cols[:11]
    quant_cols = cols[11:26] + ['FisherScore', 'EBscore', 'PoN-Ref', 'PoN-Alt']
    if config['EBFilter']['full_pon_output']:
        quant_cols.append('A|a|G|g|C|c|T|t|I|i|D|d')
    clin_cols=['cosmic70_ID', 'cosmic70_freq', 'cosmic70_type', 'cosmic70_score', 'cosmic90_MutID', 'cosmic90_type', 'cosmic90_score', 'CLNALLELEID', 'CLNDN', 'CLNSIG', 'clinvar_score', 'icgc29_ID', 'icgc29_freq', 'Clin_score']
    pop_col = cols[29:32]
    # the added extracted and score columns make up 8 columns:
    # 1-3:  'cosmic70_ID', 'cosmic70_freq', 'cosmic70_type' 
    # 4:    'icgc29_freq'
    # 5-7:  'clinvar_score', 'cosmic70_score', 'cosmic90_score'
    # 8:    'Clin_score'
    pred_col = cols[40:-9]
    new_cols = start_cols + quant_cols + clin_cols + pop_col + pred_col
    print('New columns: ', new_cols)
    return df[new_cols]


def get_clinical_scores(df):
    '''
    extract, score and realign clinical columns
    '''

    # INFERRED COLUMNS
    print('Add Cosmic70 derived columns')
    df = addCosmic70(df)
    print('Add ICGC derived columns')
    df = addICGC(df)

    print('Derive Clinical Scores from Cosmic and Clinvar')   
    # SCALAR SCORES FROM CLINICAL DBS
    df['clinvar_score'] = df.apply(get_CLINVARscore, axis=1)
    df['cosmic70_score'] = df.apply(cosmic70score, axis=1)
    df['cosmic90_score'] = df['cosmic90_type'].str.extractall(pattern).apply(cosmic90_score, axis=1).reset_index().drop(columns='match').groupby('level_0').sum()
    df['cosmic90_score'] = df['cosmic90_score'].fillna(0)

    # GET COMBINED CLIN_SCORE
    df['Clin_score'] = 0
    print('      Combining clinical scores into ClinScore')
    for col in clinscore_cols:
        print("            ", col)
        df.loc[df[col] != ".", 'Clin_score'] += ClinScore[col] * df[df[col] != "."][col]

    # RESORT THE COLUMNS
    df = resort_cols(df)
    print(list(df.columns))
    return df


def get_gene_lists(df):
    '''
    add gene info from attached lists
    '''

    ### candidate list #########
    list_cols = []
    if args.candidate_list:
        candidate_file = os.path.join(args.static_path, args.candidate_list)
        print(f'Using candidate list {candidate_file}')
        candidate_list = list(pd.read_csv(candidate_file, header=None)[0])
        df['isCandidate'] = df['Gene'].isin(candidate_list).astype(int)
        list_cols.append('isCandidate')
    if args.driver_list:
        driver_file = os.path.join(args.static_path, args.driver_list)
        print(f'Using candidate list {driver_file}')
        driver_list = list(pd.read_csv(driver_file, header=None)[0])
        df['isDriver'] = df['Gene'].isin(driver_list).astype(int)
        list_cols.append('isDriver')

    # resort columns
    cols = list(df.columns)
    start_cols = cols[:11]
    # the last two columns are the newly addded columns isDriver and isCandidate
    # have to be omitted
    rest_cols = cols[11:-2]
    new_cols = start_cols + list_cols + rest_cols
    print(new_cols)
    return df[new_cols]


print(list(anno_df.columns))
clin_df = get_clinical_scores(anno_df)

candidate_df = get_gene_lists(clin_df)
print(list(candidate_df.columns))

# this is raw unfiltered data, only informative columns were added
candidate_df.to_csv(o, sep='\t', index=False)
print(f"Writing mutation list with added columns to {o}.")


#####################################################################
############### FUNC FILTER #########################################


def filter_unwanted(df, keep_syn=False):
    '''
    basic cutoff based on gene function
    '''

    exon_func = df['ExonicFunc'] != "unknown" if keep_syn else ~df['ExonicFunc'].isin(["unknown", "synonymous SNV"])#  & (df['ExonicFunc'].notna())   
    aa_change = True # (df['AAChange'] != "UNKNOWN") & df['AAChange'].notna()  # keep for splicing
    function = ~df['Func'].isin(["downstream","intergenic","intronic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "upstream", "upstream;downstream", "UTR3", "UTR5", "UTR5;UTR3"])
    # somatic = df['somatic_status'] != 'Germline'
    return df[exon_func & aa_change & function]


basic_df = filter_unwanted(candidate_df, keep_syn=keep_syn)
basic_file = f"{output_base}.basic.csv"


basic_df.to_csv(basic_file, sep='\t', index=False)
print(f"Writing basic filtered list to {basic_file}.")

#####################################################################
# ############## BASIC FILTER #########################################

filter_setting = {
    'strict': {
        'variantT': 10,
        'Tdepth': 40,
        'EBscore': 1.4,
        'FisherScore': 20,
        'TVAF': 0.07,
        'NVAF': 0.02,
        'P-value': 2,
        'Clin_score': 100
    },
    'moderate': {
        'variantT': 3,
        'Tdepth': 30,
        'EBscore': 5,
        'FisherScore': 25,
        'TVAF': 0.05,
        'NVAF': 0.25,
        'P-value': 0.8,
        'Clin_score': 50  # value for the rescue
    },
    'loose': {
        'variantT': 10,
        'Tdepth': 20,
        'EBscore': 1,
        'FisherScore': 40,
        'TVAF': 0.01,
        'NVAF': 0.3,
        'Clin_score': 20
    },
    'Gnomad': 0.01,
    'esp6500siv2': 0.1
}


def basic_filter(df, stringency='loose'):

    # get thresholds
    thresh = filter_setting[stringency]
    tumor_depth = (df['TR2'] > thresh['variantT']) & (df['Tdepth'] > thresh['Tdepth'])
    # EBFilter
    eb = df['EBscore'] >= thresh['EBscore']

    # Strand Ratio (as FisherScore and simple)
    # strand = (df['FisherScore'] <= thresh['FisherScore']) & (df['TR2+'] != df['TR2']) & (df['TR2+'] != 0)
    strand = (df['TR2+'] != df['TR2']) & (df['TR2+'] != 0)
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])
    # VAF = (df['TVAF'] >= thresh['TVAF'])
    # total_score = df['Total_score'] >= thresh['Total_score']
    return df[tumor_depth & eb & strand & VAF]


# from Kenichi Data
# misRate_tumor > 0.05 or CoSMIC OCCURRENCE
# depth_tumor > 30 ?
# variantNum_tumor > 3
# EBscore >=4 or CoSMIC OCCURRENCE


def damm_filter(data, output):
    '''
    creates filtered output using the daniel filtering
    input: pd.dataframe of unfiltered annovar output
    output:
    - filtered/sample_tumor-normal_daniel.csv
    '''
    ############## FILTERS ##########################################

    filter_stats = {}
    filter_stats['initial'] = len(data.index)

    # filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
    exonic = filter_unwanted(data)
    filter_stats.update({'exonic':len(exonic.index)})

    exonic = get_FS_col(exonic)

    # take out reads below basic thresholds
    loose = basic_filter(exonic, 'loose')
    filter_stats.update({'loose_filter':len(loose.index)})

    moderate = basic_filter(exonic, 'moderate')
    filter_stats.update({'moderate_filter':len(moderate.index)})

    strict = basic_filter(exonic, 'strict')
    filter_stats.update({'strict_filter':len(strict.index)})

    # ########### DataBase filtering #################################
    db_loose = DB_filter(loose)
    filter_stats.update({'loose_DB_filtered':len(db_loose.index)})

    db_moderate = DB_filter(moderate)
    filter_stats.update({'moderate_DB_filtered':len(db_moderate.index)})

    db_strict = DB_filter(strict)
    filter_stats.update({'strict_DB_filtered':len(db_strict.index)})


# write selected output to files
    exonic.to_csv(f"{output}.csv", sep='\t', na_rep=".", index=False)
    moderate.to_csv(f"{output}_moderate.csv", sep='\t', na_rep=".", index=False)
    strict.to_csv(f"{output}_strict.csv", sep='\t', na_rep=".", index=False)
    db_moderate.to_csv(f"{output}_db_moderate.csv", sep='\t', na_rep=".", index=False)
    db_strict.to_csv(f"{output}_db_strict.csv", sep='\t', na_rep=".", index=False)

    return filter_stats

# print(damm_filter(df, output_base))
