#!/usr/bin/env python

import os
import argparse
import pandas as pd
import openpyxl
import primer3
from functools import partial
import math

# set_up the parser for input
parser = argparse.ArgumentParser('filters annovar output with custom criteria')
parser.add_argument('-static_path', type=str, default='', help='path to the static folder containing the gene lists')
parser.add_argument('-ref_split', type=str, default='', help='path to the split genome file used by primer3')
parser.add_argument('-keep_syn', type=str, default=False, help='True if you want to keep synonymous mutations')
parser.add_argument('-threads', type=int, default=1, help='threads for assuring memory')
parser.add_argument('input', type=str, help='input from annovar')
parser.add_argument('output', type=str, help='output file to filtered/..')


# read arguments
args = parser.parse_args()
keep_syn = args.keep_syn.lower() in ['true', 'yes', 't']
i, o = args.input, args.output
threads = args.threads
output_base = os.path.splitext(o)[0]

# genome for primer3
genome_split = os.path.join(args.static_path, args.ref_split)

print(f'Started editing and basic filtering for {i}.')
anno_df = pd.read_csv(i, sep='\t')

print(f"keep_syn= {keep_syn}")


#  ############## BASIC FILTER ####################################
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
    'PoN-Alt-NonZeros': 4,
    'FisherScore': 50,
    'TVAF': 0.01,
    'NVAF': 0.3,
}


def filter1(df):

    # get thresholds
    thresh = filter1_setting
    tumor_depth = (df['TR2'] > thresh['variantT']) & (df['Tdepth'] > thresh['Tdepth'])

    # EBFilter
    if thresh['EBscore']:
        eb = df['EBscore'] >= thresh['EBscore']
    else:
        eb = True

    # other PanelOfNormal metrices
    pon = (df['PoN-Ratio'] < thresh['PoN-Ratio']) | (df['PoN-Alt-NonZeros'] <= thresh['PoN-Alt-NonZeros'])

    # Strand Bias (as FS)
    no_strand_bias = df['FisherScore'] < thresh['FisherScore']
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF']) & (df['TVAF'] > thresh['NVAF'])

    return df[tumor_depth & (eb | pon) & no_strand_bias & VAF]


# from Kenichi Data
# misRate_tumor > 0.05 or CoSMIC OCCURRENCE
# depth_tumor > 30 ?
# variantNum_tumor > 3
# EBscore >=4 or CoSMIC OCCURRENCE

filter1_file = f"{output_base}.filter1.csv"
filter1_df = filter1(basic_df)

########################################################################
# ##### Primer3 ########################################################

PCR_config = {
    'seq_len': 500,
    'mut_pad': 5,
    'prod_size_min': 140,
    'prod_size_max': 200
}

Primer3_config = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    }


def file2str(file):
    '''
    returns a string from a text file
    '''

    with open(file, 'r') as file:
        return file.read().upper().replace('\n', '')


def get_primer_df(chrom_seq, config, row):
    '''
    returns the best primer pair for a given position
    return value is [fwd_seq, fwd_tmp, rev_seq, rev_tmp, prod_size]
    active chromosome sequence is global variable chrom
    '''

    # load sequence
    pos = row['Start']
    half_seq = int(config['seq_len'] / 2)
    seq_start = pos - half_seq
    seq_end = pos + half_seq
    seq = chrom_seq[seq_start:seq_end]
    pad = int(config['mut_pad'] / 2)
    half_size = int(config['prod_size_min'] / 2)

    # calculate the target_range as offSet from center (half)
    offSet = half_size - 20 - pad
    target_start = half_seq - offSet
    target = [target_start, offSet * 2]
    setting = {
        'SEQUENCE_ID': 'asdf',
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_TARGET': target
    }
    primers = primer3.bindings.designPrimers(setting, config)

    # return '--' if nothing was found
    if primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
        row['fwd_seq'] = row['fwd_tmp'] = row['rev_seq'] = row['rev_tmp'] = row['prod_size'] = '--'
        return row

    prod_len = primers['PRIMER_RIGHT_0'][0] - primers['PRIMER_LEFT_0'][0]
    row['PrimerFwd'] = primers['PRIMER_LEFT_0_SEQUENCE']
    row['PrimerFwd_Tmp'] = int(primers['PRIMER_LEFT_0_TM'] * 100) / 100
    row['PrimerRev'] = primers['PRIMER_RIGHT_0_SEQUENCE']
    row['PrimerRev_Tmp'] = int(primers['PRIMER_RIGHT_0_TM'] * 100) / 100
    row['ProductSize'] = prod_len
    row['MutLocation'] = half_seq - primers['PRIMER_LEFT_0'][0]

    return row


def run_primer3(mut_df, genome_split_folder='', pcr_config=PCR_config, primer3_config=Primer3_config):

    # apply pcr size to primer3_config
    primer3_config['PRIMER_PRODUCT_SIZE_RANGE'] = [pcr_config['prod_size_min'], pcr_config['prod_size_max']]
    primer3_config.update(pcr_config)

    mut_df.loc[:, 'Chr'] = mut_df['Chr'].astype('str')
    org_cols = list(mut_df.columns)
    df_list = []
    # cycle through (formatted) chromosomes
    # + load chromosome sequence
    # + create primer_df for mutations on that chromosome
    # + concat all mutations
    for chrom in mut_df['Chr'].unique():
        chrom_seq = file2str(f'{genome_split_folder}/{chrom}.fa')
        chr_df = mut_df.query('Chr == @chrom')
        primer_df = chr_df.apply(partial(get_primer_df, chrom_seq, primer3_config), axis=1)
        df_list.append(primer_df)
    primer_df = pd.concat(df_list, sort=True)
    clinscore_index = list(mut_df.columns).index('Clin_score') + 1
    primer_df = primer_df[org_cols[:clinscore_index] + ['PrimerFwd', 'PrimerFwd_Tmp', 'PrimerRev', 'PrimerRev_Tmp', 'ProductSize', 'MutLocation'] + org_cols[clinscore_index:]]
    return primer_df


filter1_p3 = run_primer3(filter1_df, genome_split_folder=genome_split, pcr_config=PCR_config, primer3_config=Primer3_config)

filter1_p3.to_csv(filter1_file, sep='\t', index=False)
len_filter1 = len(filter1_p3.index)
print(f"Writing filter1 list ({len_filter1} mutations) to {filter1_file}.")


filter2_setting = {
    'strict': {
        'variantT': 10,
        'Tdepth': 60,
        'EBscore': 5,
        'PoN-Ratio': 0,
        'PoN-Alt-NonZeros': 1,
        'FisherScore': 20,
        'strand_polarity':1, # filters out TR2 == X and (TR2+ <= 1 | TR2+ >= X-1)
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
        'PoN-Alt-NonZeros': 2,
        'FisherScore': 30,
        'strand_polarity': 0,
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
        'PoN-Alt-NonZeros': 3,
        'FisherScore': 40,
        'strand_polarity': None,
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
    if thresh['EBscore']:
        eb = df['EBscore'] >= thresh['EBscore']
    else:
        eb = True

    pon = (df['PoN-Ratio'] < thresh['PoN-Ratio']) | (df['PoN-Alt-NonZeros'] < thresh['PoN-Alt-NonZeros'])

    # minimum TVAF if not candidate
    is_candidate = (df['isCandidate'] == 1) | (df['isDriver'] == 1)
    no_noise = is_candidate | (df['TVAF'] > 0.05)

    # Strand Ratio (as FisherScore and simple)
    no_strand_bias = df['FisherScore'] <= thresh['FisherScore']

    # Strand Polarity (filters out very uneven strand distribution of alt calls)
    if thresh['strand_polarity']:
        pol = thresh['strand_polarity']
        no_strand_polarity = (df['TR2+'] <= df['TR2'] - pol) & (df['TR2+'] >= pol)
    else:
        no_strand_polarity = True

    strandedness = no_strand_bias & no_strand_polarity

    # VAF is simple
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])

    # Clin_score is used for rescue of all mutations
    clin_score = df['Clin_score'] >= thresh['Clin_score']

    # apply filters to dataframe
    df = df[(tumor_depth & (pon | eb) & strandedness & VAF & no_noise) | clin_score].sort_values(['TVAF'], ascending=False)
    list_len = len(df.index)
    return df, list_len


print(f"Applying filter2 in 3 stringencies [loose, moderate, strict]")
df_filter2_loose, len_loose = filter2(filter1_p3, stringency='loose')
df_filter2_mod, len_mod = filter2(filter1_p3, stringency='moderate')
df_filter2_strict, len_strict = filter2(filter1_p3, stringency='strict')

print('loose:', len_loose)
print('moderate:', len_mod)
print('strict:', len_strict)

df_filter2_loose.to_csv(f"{output_base}.loose.csv", sep='\t', index=False)
df_filter2_mod.to_csv(f"{output_base}.moderate.csv", sep='\t', index=False)
df_filter2_strict.to_csv(f"{output_base}.strict.csv", sep='\t', index=False)
print(f"Writing filter2 lists to {output_base}.<stringency>.csv")


excel_file = f"{output_base}.filter.xlsx"
with pd.ExcelWriter(excel_file) as writer:
    filter1_p3.to_excel(writer, sheet_name=f'filter1 <{len_filter1}>', index=False)
    df_filter2_loose.to_excel(writer, sheet_name=f'loose <{len_loose}>', index=False)
    df_filter2_mod.to_excel(writer, sheet_name=f'moderate <{len_mod}>', index=False)
    df_filter2_strict.to_excel(writer, sheet_name=f'strict <{len_strict}>', index=False)


print(f"Writing combined filters to excel file {excel_file}.")

