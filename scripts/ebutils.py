from ebcore import matrix2EBscore, AB2EBscore, matrix2AB
from multiprocessing import Pool
from functools import partial
from datetime import datetime as dt
import os
import numpy as np
import math
import pandas as pd

ansii_colors = {
          'magenta': '[1;35;2m',
          'green': '[1;9;2m',
          'red': '[1;31;1m',
          'cyan': '[1;36;1m',
          'gray': '[1;30;1m',
          'black': '[0m'
          }

colors = {
        'process': ansii_colors['green'],
        'time': ansii_colors['magenta'],
        'normal': ansii_colors['gray'],
        'warning': ansii_colors['red'],
        'success': ansii_colors['cyan']
        }


def show_output(text, color='normal', multi=False, time=True):
    '''
    get colored output to the terminal
    '''
    time = f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ''
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ''
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text)


def show_command(command, list=False, multi=True):
    '''
    prints the command line if debugging is active
    '''

    proc = f"\033[92mProcess {os.getpid()}\033[0m : " if multi else ""
    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(proc + command)
    return


# # check if sample is included in PoN ######
def checkPon4sample(row, bam_base):
    '''
    per row of Pon_list checks whether sample name occurs in Pon_list
    if yes, returns the row number (zero-based)
    '''
    if bam_base in row.iloc[0]:
        print(row[0], row.name)
        return row.name
    return None


def get_sample_pos(pon_list, bam_file):
    '''
    returns the zero-based position of corresponding normal sample in Pon_list
    '''
    # load pon_list as df
    pon_df = pd.read_csv(pon_list, header=None)
    # get sample base_name
    bam_base = os.path.basename(bam_file).split('.')[0].split("_")[0]

    sample_pos = pon_df.apply(checkPon4sample, axis=1, bam_base=bam_base).sum()
    return int(sample_pos + 1)
# ######### ADD BASE INFO ##############################################


def get_pon_bases(matrix_df, remove_sample=True):
    '''
    returns from eb-matrix file the concatenated Pon coverage for pos and neg strand
    this is important output for mutation QC
    imput cols:
        depthP
        depthN
        misP
        misN
    '''
    if remove_sample:
        # remove sample depths from the columns
        for col in ['depthP', 'misP', 'depthN', 'misN']:
            matrix_df[col] = matrix_df[col].str.replace(r"^[0-9]+\|", "")

    # concate the respective columns
    matrix_df['PoN-Ref'] = matrix_df['depthP'].str.cat(matrix_df['depthN'], sep="-")
    matrix_df['PoN-Alt'] = matrix_df['misP'].str.cat(matrix_df['misN'], sep="-")
    return matrix_df


# ################# EB from matrix ###########################################
def compute_matrix2EB(fit_pen, df):
    '''
    per df row, computes the EBscore from full depth matrix
    first row: target depth
    next rows: pon depth
    '''

    show_output(f"Computing EBscore for {len(df.index)} lines", multi=True)
    df['EBscore'] = df.apply(partial(matrix2EBscore, fit_pen), axis=1)
    show_output("Finished!", multi=True, time=True)
    return df


def compute_matrix2EB_multi(df, pen, threads):
    print('Pooling threads: ', threads)
    eb_pool = Pool(threads)
    print('Pooled')
    # minimal length of 1000 lines
    split_factor = min(math.ceil(len(df.index) / 1000), threads)
    split = np.array_split(df, split_factor)
    dfs = eb_pool.map(partial(compute_matrix2EB, pen), split)

    # out_df contains EB_score
    out_df = pd.concat(dfs).sort_values(['Chr', 'Start'])
    return out_df

# ################# EB from AB ###########################################


def compute_AB2EB(df):
    '''
    per row of df, takes a target depth-ponAB matrix and computes the EBscore 
    '''

    show_output(f"Computing EBscore for {len(df.index)} lines", multi=True)
    df['EBscore'] = df.apply(AB2EBscore, axis=1)
    show_output("Finished!", multi=True)
    return df


def compute_AB2EB_multi(df, threads):

    eb_pool = Pool(threads)
    # minimal length of 2000 lines
    split_factor = min(math.ceil(len(df.index) / 2000), threads)
    df_split = np.array_split(df, split_factor)
    dfs = eb_pool.map(compute_AB2EB, df_split)
    # out_df contains EB_score
    out_df = pd.concat(dfs).sort_values(['Chr', 'Start'])
    return out_df


# ##################### EB-CACHE ###########################

def computeEBcache(mat_df, pen):
    show_output(f"Computing EBcache for {len(mat_df.index)} lines", time=True, multi=True, color='process')
    cache_df = matrix2AB(mat_df, pen)
    show_output(f"Finished!", time=True, multi=True, color='success')
    return cache_df


def matrix2AB_multi(matrix_file, output, pen, threads):
    matrix_df = pd.read_csv(matrix_file, sep='\t', compression='gzip', index_col=False)
    cache_pool = Pool(threads)
    matrix_split = np.array_split(matrix_df, threads)
    cache_dfs = cache_pool.map(computeEBcache, matrix_split)
    cache_df = pd.concat(cache_dfs)
    cache_df.to_csv(output, compression='gzip', sep='\t', index=False)
    return output