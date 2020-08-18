from ebcore import matrix2EBscore, AB2EBscore, matrix2AB
from multiprocessing import Pool
from functools import partial
from itertools import repeat

import os
import numpy as np
import math
import pandas as pd
from script_utils import show_output


# # check if sample is included in PoN ######
def checkPon4sample(row, bam_base):
    '''
    per row of Pon_list checks whether sample name occurs in Pon_list
    if yes, returns the row number (zero-based)
    '''

    # the split is necessary for looking only in file and not in paths
    if bam_base in row.iloc[0].split('/')[-1]:
        print(row[0], row.name)
        return int(row.name) + 1
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
    return int(sample_pos)
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
    matrix_df['PoN-Ref'] = matrix_df['depthP'].str.cat(
        matrix_df['depthN'], sep="-")
    matrix_df['PoN-Alt'] = matrix_df['misP'].str.cat(
        matrix_df['misN'], sep="-")
    return matrix_df


# ################# EB from matrix ###########################################
def compute_matrix2EB(df, fit_pen):
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
    '''
    split --> Pool --> compute --> concate
    '''

    eb_pool = Pool(threads)

    # minimal length of 200 lines
    split_factor = min(math.ceil(len(df.index) / 200), threads)
    split = np.array_split(df, split_factor)
    dfs = eb_pool.starmap(compute_matrix2EB, zip(split, repeat(pen)))

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
    show_output(
        f"Computing EBcache for {len(mat_df.index)} lines", time=True, multi=True, color='process')
    cache_df = matrix2AB(mat_df, pen)
    show_output(f"Finished!", time=True, multi=True, color='success')
    return cache_df


def matrix2AB_multi(matrix_file, output, pen, threads):
    matrix_df = pd.read_csv(matrix_file, sep='\t',
                            compression='gzip', index_col=False)
    cache_pool = Pool(threads)
    matrix_split = np.array_split(matrix_df, threads)
    # in order to map both mat_df and pen into computeEBcache
    # pool.starmap has to be used with an iterable of argument tuples.
    # argument tuples can infused with the zip(a, repeat(b))
    # pattern for the constant pen variable
    cache_dfs = cache_pool.starmap(
        computeEBcache,
        zip(matrix_split, repeat(pen))
    )
    cache_df = pd.concat(cache_dfs)
    cache_df.to_csv(output, compression='gzip', sep='\t', index=False)
    return output
