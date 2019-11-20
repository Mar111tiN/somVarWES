#!/usr/bin/env python

import argparse
import pandas as pd

# set_up the parser for input
parser = argparse.ArgumentParser('Merges input from annovar to merged csv')
parser.add_argument('-indel', '-i', type=str, help='indel output from annovar')
parser.add_argument('-snp', '-s', type=str, help='snp output from annovar')
parser.add_argument('-out', '-o', type=str, help='name of merged csv')


# read arguments
args = parser.parse_args()
i, s, o = args.indel, args.snp, args.out

# read in data, add respective column and add to data


def convert_vaf(vaf):
    '''
    converts percent VAF to frequency VAF
    '''
    try:
        return (float(vaf.split('%')[0]) / 100)
    except:
        return vaf


data = []
for mut_type in [(i, 'indel'), (s, 'snp')]:
    mut_df = pd.read_csv(mut_type[0], sep='\t', dtype={'Chr': str}, header=0, converters={'NVAF': convert_vaf, 'TVAF': convert_vaf}, na_values='.')
    rows = len(mut_df.index)
    # create the Series of correct length containing the mut_type as string 
    mut_type = pd.Series(mut_type[1], index=range(rows), name='mut_type')
    mut_labeled = pd.concat([mut_df, mut_type], axis=1)
    data.append(mut_labeled)


# merge the two files and write to data
anno_merge = pd.concat(data, axis=0).sort_values(['Chr', 'Start'], ascending=True)
anno_merge['Tdepth'] = anno_merge['TR1'] + anno_merge['TR2']
anno_merge['Ndepth'] = anno_merge['NR1'] + anno_merge['NR2']
# resort the columns and omit Tgenotype and Ngenotype
col = anno_merge.columns
first_cols = list(col[:10])
varscan_cols = ['mut_type', 'somatic_status', 'variantP', 'somaticP', 'NVAF', 'Ndepth', 'NR1', 'NR1+', 'NR2', 'NR2+', 'TVAF', 'Tdepth', 'TR1', 'TR1+', 'TR2', 'TR2+']
anno_cols = list(col[10:-22])
merged_cols = first_cols + varscan_cols + anno_cols
anno_merge = anno_merge[merged_cols]
anno_merge.to_csv(o, index=False, sep='\t')