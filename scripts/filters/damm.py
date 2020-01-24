#!/usr/bin/env python

import os
import argparse
import pandas as pd
import openpyxl

# set_up the parser for input
parser = argparse.ArgumentParser('filters annovar output with custom criteria')
parser.add_argument('-static_path', type=str, default='', help='path to the static folder containing the gene lists')
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

print(f'Started editing and basic filtering for {i}.')
anno_df = pd.read_csv(i, sep='\t')


# ############## BASIC FILTER ####################################
def filter_basic(df, keep_syn=False):
    '''
    basic cutoff based on gene function
    '''

    exon_func = df['ExonicFunc'] != "unknown" if keep_syn else ~df['ExonicFunc'].isin(["unknown", "synonymous SNV"])  #  & (df['ExonicFunc'].notna())   
    aa_change = True  # (df['AAChange'] != "UNKNOWN") & df['AAChange'].notna()  # keep for splicing
    function = ~df['Func'].isin(["downstream", "intergenic", "intronic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "upstream", "upstream;downstream", "UTR3", "UTR5", "UTR5;UTR3"])
    # somatic = df['somatic_status'] != 'Germline'
    return df[exon_func & aa_change & function]


basic_df = filter_basic(anno_df, keep_syn=keep_syn)
basic_file = f"{output_base}.csv"


basic_df.to_csv(basic_file, sep='\t', index=False)
print(f"Writing basic filtered list to {basic_file}.")

# ############### FILTER1 ########################################

filter1_setting = {
    'variantT': 2,
    'Tdepth': 20,
    'EBscore': 1,
    'PoN-Ratio': 0.001,
    'PoN-Alt-Zeros': 46,
    'FisherScore': 50,
    'TVAF': 0.01,
    'NVAF': 0.3,
}


def filter1(df):

    # get thresholds
    thresh = filter1_setting
    tumor_depth = (df['TR2'] > thresh['variantT']) & (df['Tdepth'] > thresh['Tdepth'])

    # EBFilter
    # eb = df['EBscore'] >= thresh['EBscore']
    eb = (df['PoN-Ratio'] < thresh['PoN-Ratio']) | (df['PoN-Alt-Zeros'] > thresh['PoN-Alt-Zeros'])

    # Strand Ratio (as FisherScore and simple)
    strand = (df['TR2+'] < df['TR2']) & (df['TR2+'] > 0)
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])
    return df[tumor_depth & eb & strand & VAF]


# from Kenichi Data
# misRate_tumor > 0.05 or CoSMIC OCCURRENCE
# depth_tumor > 30 ?
# variantNum_tumor > 3
# EBscore >=4 or CoSMIC OCCURRENCE

filter1_file = f"{output_base}.filter1.csv"
filter1_df = filter1(basic_df)
filter1_df.to_csv(filter1_file, sep='\t', index=False)
print(f"Writing filter1 list to {filter1_file}.")


filter2_setting = {
    'strict': {
        'variantT': 10,
        'Tdepth': 60,
        'EBscore': 5,
        'PoN-Ratio': 0,
        'PoN-Alt-Zeros': 49,
        'FisherScore': 20,
        'TVAF': 0.03,
        'NVAF': 0.25,
        'P-value': 2,
        'Clin_score': 500
    },
    'moderate': {
        'variantT': 7,
        'Tdepth': 50,
        'EBscore': 4,
        'PoN-Ratio': 0.0001,
        'PoN-Alt-Zeros': 48,
        'FisherScore': 30,
        'TVAF': 0.02,
        'NVAF': 0.35,
        'P-value': 0.8,
        'Clin_score': 100
    },
    'loose': {
        'variantT': 5,
        'Tdepth': 40,
        'EBscore': 2,
        'PoN-Ratio': 0.0005,
        'PoN-Alt-Zeros': 47,
        'FisherScore': 40,
        'TVAF': 0.01,
        'NVAF': 0.4,
        'Clin_score': 50
    },
    'Gnomad': 0.01,
    'esp6500siv2': 0.1
}


def filter2(df, stringency='moderate'):
    '''
    creates filtered output using the daniel filtering
    input: pd.dataframe of unfiltered annovar output
    output:
    - filtered/sample_tumor-normal_daniel.csv
    '''
    # get thresholds
    thresh = filter2_setting[stringency]
    tumor_depth = (df['TR2'] > thresh['variantT']) & (df['Tdepth'] > thresh['Tdepth'])
    # eb = df['EBscore'] >= thresh['EBscore']
    eb = (df['PoN-Ratio'] < thresh['PoN-Ratio']) | (df['PoN-Alt-Zeros'] > thresh['PoN-Alt-Zeros'])
    # minimum TVAF if not candidate
    is_candidate = (df['isCandidate'] == 1) | (df['isDriver'] == 1)
    no_noise = is_candidate | (df['TVAF'] > 0.05)

    # Strand Ratio (as FisherScore and simple)
    strand = df['FisherScore'] <= thresh['FisherScore']
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])
    clin_score = df['Clin_score'] >= thresh['Clin_score']
    df = df[(tumor_depth & eb & strand & VAF & no_noise) | clin_score].sort_values(['TVAF'], ascending=False)
    list_len = len(df.index)
    return df, list_len


df_filter2_loose, len_loose = filter2(filter1_df, stringency='loose')
df_filter2_mod, len_mod = filter2(filter1_df, stringency='moderate')
df_filter2_strict, len_strict = filter2(filter1_df, stringency='strict')

print('loose:', len_loose)
print('moderate:', len_mod)
print('strict:', len_strict)

df_filter2_loose.to_csv(f"{output_base}.loose.csv", sep='\t', index=False)
df_filter2_mod.to_csv(f"{output_base}.moderate.csv", sep='\t', index=False)
df_filter2_strict.to_csv(f"{output_base}.strict.csv", sep='\t', index=False)

with pd.ExcelWriter(f"{output_base}.filter2.xlsx") as writer:
    df_filter2_loose.to_excel(writer, sheet_name=f'loose <{len_loose}>', index=False)
    df_filter2_mod.to_excel(writer, sheet_name=f'moderate <{len_mod}>', index=False)
    df_filter2_strict.to_excel(writer, sheet_name=f'strict <{len_strict}>', index=False)
