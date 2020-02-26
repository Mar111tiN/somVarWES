import primer3
import re
import math
import pandas as pd

import os
from functools import partial


############# SNAKEMAKE ##################
w = snakemake.wildcards
config = snakemake.config
p_config = config['primer3']

i = snakemake.input[0]
o = snakemake.output[0]
threads = snakemake.threads
genome_split = snakemake.params.genome_split
# ##### Primer3 ########################################################

PCR_config = {
    'seq_len': 500,
    'mut_pad': 5,
    'prod_size_min': p_config['min_size'],
    'prod_size_max': p_config['max_size']
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


def get_chrom(chrom, chroms_folder='chroms'):
    '''
    convenience function returning the chromosome sequence without
    the header sequences
    a folder containing the split chromosomes has to be provided
    '''

    if 'CHR' in str(chrom).upper():
        chrom = chrom.upper().replace('CHR', '')
    chrom_file = f"{chroms_folder}/chr{chrom}.fa"
    chrom_seq = file2str(chrom_file)
    chrom_seq = re.sub(r'^>CHR[0-9XY]+', '', chrom_seq)
    chrom = {
        "name": f"chr{chrom}",
        "sequence": chrom_seq,
        "length": len(chrom_seq)
    }
    return chrom 


def mut2insert(mut={}, seq={}, return_none=False):
    '''
    takes a mutation dictionary of shape:
    {
        'Chrom': 'chr4,
        'Start': 13331,
        'End': 13331,
        'Ref': '-',
        'Alt': 'A'
        }
    and a seq dictionary of shape:
    {
        'Chrom': 'chr4,
        'Start': 13300,
        'End': 123123,
        'seq': 'ATTTCTCCCACTCCCACA',
        }
    and returns the seq with the mutation inserted into the sequence
    if mutation location is out of bounds of sequence, only the sequence is returned without editing
    '''

    if mut['Chr'] != seq['Chr']:
        return None if return_none else seq['seq']
    if (mut['Start'] < seq['Start']) or (mut['End'] > seq['End']):
        return None if return_none else seq['seq']
    bases = ['A', 'C', 'G', 'T']
    # case SNP
    start = mut['Start'] - seq['Start']
    end = mut['End'] - seq['Start']

    if mut['Ref'] == "-":
        upstream = seq['seq'][:start]
        downstream = seq['seq'][end:]
        mutation = f"<<+{mut['Alt']}>>"
    else:
        upstream = seq['seq'][:start]
        downstream = seq['seq'][end + 1:]
        if mut['Alt'] == "-":
            mutation = f">∆{mut['Ref']}<"
        else:
            # SNP
            mutation = f"({mut['Ref']}>{mut['Alt']})"
    return f"{upstream}{mutation}{downstream}"


def edit_seq(row):
    # convert mutation to dict
    def mut2dict(row):
        # mutation to dict
        chrom = row['Chr']
        start = row['Start']
        end = row['End']
        mut_dict = {
            'Chr': row['Chr'],
            'Start': row['Start'],
            'End': row['End'],
            'Ref': row['Ref'],
            'Alt': row['Alt']
        }
        return mut_dict
    mut_dict = mut2dict(row)


    insert_dict = {
        'Chr': row['Chr'],
        'Start': row['InsertStart'],
        'End': row['InsertEnd'],
        'seq': row['InsertSeq']
    }

    edited_seq = mut2insert(mut=mut_dict, seq=insert_dict, return_none=True)
    return edited_seq


def get_primer_df(chrom, config, row):
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
    seq = chrom['sequence'][seq_start:seq_end]
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
    
    ## get chrom coords
    amp_start = seq_start + primers['PRIMER_LEFT_0'][0] + 1
    amp_end = seq_start + primers['PRIMER_RIGHT_0'][0] + 1
    
    row['AmpliconRange'] = f"{chrom['name']}:{amp_start}-{amp_end}"
    
    insert_start = amp_start + primers['PRIMER_LEFT_0'][1]
    insert_end = amp_end - primers['PRIMER_RIGHT_0'][1]
    row['InsertRange'] = f"{chrom['name']}:{insert_start}-{insert_end}"
    row['InsertSize'] = insert_end - insert_start
    insert_seq = chrom['sequence'][insert_start -1:insert_end]
    row['InsertSeq'] = insert_seq
    
    mut_dict = {
        'Chr': row['Chr'],
        'Start': row['Start'],
        'End': row['End'],
        'Ref': row['Ref'],
        'Alt': row['Alt']
    }
    
    insert_dict = {
        'Chr': row['Chr'],
        'Start': insert_start,
        'End': insert_end,
        'seq': insert_seq
    }
    
    
    row['InsertSeq'] = mut2insert(mut=mut_dict, seq=insert_dict, return_none=True)
    row['offsetL'] = row['Start'] - insert_start
    row['offsetR'] = insert_end - row['End']
    
    
    
    row['fwd-Primer'] = primers['PRIMER_LEFT_0_SEQUENCE']
    row['rev_Primer'] = primers['PRIMER_RIGHT_0_SEQUENCE']
    row['AmpliconSize'] = sum(primers['PRIMER_RIGHT_0']) - primers['PRIMER_LEFT_0'][0]
    row['Status'] = 'not established'
    row['Note'] = f"(fwd={int(primers['PRIMER_LEFT_0_TM'] * 10) / 10}|rev={int(primers['PRIMER_RIGHT_0_TM'] * 10) / 10})"
    return row


def run_primer3(mut_df, chroms_folder='', pcr_config=PCR_config, primer3_config=Primer3_config, primer_list=''):

    # apply pcr size to primer3_config
    primer3_config['PRIMER_PRODUCT_SIZE_RANGE'] = [pcr_config['prod_size_min'], pcr_config['prod_size_max']]
    primer3_config.update(pcr_config)

    mut_df.loc[:, 'Chr'] = mut_df['Chr'].astype('str')
    mut_df = mut_df.loc[:, ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']]
    org_cols = list(mut_df.columns)
    df_list = []
    # cycle through (formatted) chromosomes
    # + load chromosome sequence
    # + create primer_df for mutations on that chromosome
    # + concat all mutations
    for chrom in mut_df['Chr'].unique():
        chromo = get_chrom(chrom, chroms_folder=chroms_folder)
        chr_df = mut_df.query('Chr == @chrom')
        primer_df = chr_df.apply(partial(get_primer_df, chromo, primer3_config), axis=1)
        df_list.append(primer_df)
    primer_df = pd.concat(df_list, sort=True)
    new_cols = [
        'fwd-Primer', 
        'rev_Primer',
        'Status',
        'Note',
        'AmpliconRange',
        'AmpliconSize',
        'InsertRange',
        'InsertSize',
        'InsertSeq',
        'offsetL',
        'offsetR'
    ]
    for col in ['AmpliconSize', 'InsertSize', 'offsetL', 'offsetR']:
        primer_df[col] = primer_df[col].fillna(0).astype(int)

    primer_df = primer_df[org_cols + new_cols]
    return primer_df
    


filter1_df = pd.read_csv(i, sep='\t')


primer_df = run_primer3(
    filter1_df,
    chroms_folder=genome_split,
    pcr_config=PCR_config,
    primer3_config=Primer3_config
)

primer_df.to_csv(o, sep='\t', index=False)
print(f"Writing primer list to {o}.")
