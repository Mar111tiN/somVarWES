import os
from os import system as shell
from os.path import getsize as filesize
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
import math
from scipy.optimize import fmin_l_bfgs_b as minimize_func
from scipy.stats import chi2
from scipy.special import gammaln
from datetime import datetime as dt


###########################################################
###################### FUNCTIONS ##########################

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

def show_output(text, color='normal', multi=False, time=False):
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
    
############### EBcache using matrix2EBinput.mawk #######################
# the matrices for beta-binomial calculation
KS_matrix = np.array([[1,0,1,1,0,1,0,0,0],[0,1,-1,0,1,-1,0,0,0]])
gamma_reduce = np.array([1,-1,-1,-1,1,1,1,-1,-1])

def bb_loglikelihood(params, count_df, is_1d):
    [a, b] = params
    ab_matrix = np.array([1,1,1,a+b,a,b,a+b,a,b])
    # convert df into matrix for np.array operations that change dims
    count_matrix = count_df.values
    # perform matrix multiplication to get inputs to log-gamma
    input_matrix = np.matmul(count_matrix,KS_matrix) + ab_matrix
    # get corresponding log-gamma values and reduce over pon-values
    if is_1d: # check whether gammatrix is 2-dim - otherwise sum aggregation over axis 0 is faulty
        gamma_matrix = gammaln(input_matrix)
    else:  
        gamma_matrix = np.sum(gammaln(input_matrix), axis=0)
    # add or subtract using gamma_reduce matrix and sum to loglikelihood (scalar)
    log_likelihood = np.sum(gamma_matrix * gamma_reduce)
    return log_likelihood
 

def fit_bb(count_df, pen):
    '''
    Obtaining maximum likelihood estimator of beta-binomial distribution
    count_df is the array of depth-mismatch (trials, success) pairs over the PoN list for either strand
    during minimization of fitting function (max for loglikelihood) penalty term is applied to constrain alpha and beta
        Ref for L-BFGS-B algorithm:
        A Limited Memory Algorithm for Bound Constrained Optimization
        R. H. Byrd, P. Lu and J. Nocedal. , (1995), 
        SIAM Journal on Scientific and Statistical Computing, 16, 5, pp. 1190-1208.
    '''

    def bb_loglikelihood_fitting(params, count_df, penalty):
        '''
        Fitting params [alpha, beta] to maximize loglikelihood
        '''

        # Here, we apply the penalty term of alpha and beta (default 0.5 is slightly arbitray...)
        result = 0.5 * math.log(sum(params)) - bb_loglikelihood(params, count_df, False) # matrix is dim2
        return result
        
    # get the respective control matrices (as dataframe) for positive and negative strands
    count_p = count_df.loc[:, ['depth_p', 'mm_p']]
    count_n = count_df.loc[:, ['depth_n', 'mm_n']]
    # minimize loglikelihood using L-BFGS-B algorithm
    ab_p = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args = (count_p, pen), approx_grad = True,
                           bounds = [(0.1, 10000000), (1, 10000000)]
                          )[0]
    ab_p = [round(param, 5) for param in ab_p]
    ab_n = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args = (count_n, pen), approx_grad = True,
                           bounds = [(0.1, 10000000), (1, 10000000)]
                          )[0]
    ab_n = [round(param, 5) for param in ab_n]
    # print(f'abP: {ab_p} - abN: {ab_n}')
    return {'p':ab_p, 'n':ab_n}


def get_count_df(row, var):
    '''
    converts the base-wise read coverage to a matrix
    '''
    
    matrix = pd.DataFrame()
    matrix['depth_p'] = np.array(row['Depth-ACGT'].split('|')).astype(int)
    matrix['mm_p'] = np.array(row[var].split('|')).astype(int)
    matrix['depth_n'] = np.array(row['Depth-acgt'].split('|')).astype(int)
    matrix['mm_n'] = np.array(row[var.lower()].split('|')).astype(int)
    if var == 'I':
        matrix['depth_p'] += np.array(row['Depth-INDEL'].split('|')).astype(int)
        matrix['depth_n'] += np.array(row['Depth-indel'].split('|')).astype(int)
    return matrix


def matrix2AB(matrix_df):
    '''
    creates the AB_df for a pileup
    '''

    AB_df = matrix_df.iloc[:, [0, 1]].copy()
    # ###################### AB FITTING ###############################################

    def get_AB(penalty, row):
        '''
        returns the AB parameters (A+a A+b A-a A-b G+a G+b....) for each pileup row
        main computational load
        '''

        bb_s = pd.Series()
        # ########## get count matrix ###########################################
        for var in list('ACTGI'):
            # get the count matrix
            count_df = get_count_df(row, var)
            # get the AB parameters for

            # <<<<<<######### DEBUG ###############
            # if row['Start'] in [
            #     14830116, 14622123, 14841205, 14733419, 14840756, 14719431, 16618785,
            #     14824039, 14830279, 14622390, 18769493, 18550359, 13905782, 14840589,
            #     14830168, 16618629, 13705217, 14830580, 18026171, 14841406, 14622699,
            #     16618619, 14719675, 13905767, 13902927, 14840574, 18026160, 14824292, 14840574, 
            #     ]:
            #     print(row, var) 
            # <<<<<<###############################


            bb_params = fit_bb(count_df, penalty)
        # dump the different parameters into bb_s
        # keys have to fit with the var_columns for the receiving AB_df
            bb_s[f'{var}+a'] = bb_params['p'][0]
            bb_s[f'{var}+b'] = bb_params['p'][1]
            bb_s[f'{var}-a'] = bb_params['n'][0]
            bb_s[f'{var}-b'] = bb_params['n'][1]
        return bb_s

    # ################# Store AB data into df ####################################
    # create the columns (A+a A+b A-a A-b G+a G+b....) for the recipient df of the pileup_df apply function
    var_columns = [f'{var}{strand}{param}' for var in list('ACTGI') for strand in ['+', '-'] for param in ['a', 'b']]

    AB_df[var_columns] = matrix_df.apply(partial(get_AB, pen), axis=1)
    return AB_df


######################### MULTITHREADING #####################################
def  computeEBcache(mat_df):
    show_output(f"Computing EBcache for {len(mat_df.index)} lines", time=True, multi=True, color='process')
    cache_df = matrix2AB(mat_df)
    show_output(f"Finished!", time=True, multi=True, color='success')
    return cache_df


def matrix2AB_multi(matrix_file, output, threads):
    matrix_df = pd.read_csv(matrix_file, sep='\t', compression='gzip', index_col=False)
    cache_pool = Pool(threads)
    matrix_split = np.array_split(matrix_df, threads)
    cache_dfs = cache_pool.map(computeEBcache, matrix_split)
    cache_df = pd.concat(cache_dfs)
    cache_df.to_csv(output, compression='gzip', sep='\t', index=False)
    return output

##############################################################################
########################## SNAKE PARAMETERS ################################s##

w = snakemake.wildcards
config = snakemake.config
threads = snakemake.threads
log = snakemake.log
cache_file = snakemake.output[0]
params = snakemake.params
pon_list = params.pon_list
EBparams = snakemake.config['EBFilter']['params']
pen = EBparams['fitting_penalty']
matrix_file = snakemake.input[0]
chrom = w.chrom
i = w.i


############### MATRIX FILE --> EBcache ##################

show_output(f"Performing AB-parameter computation for split {i} of {chrom}", color='normal', time=True)
if filesize(matrix_file) > 20:
    cache_file = matrix2AB_multi(matrix_file, cache_file, threads)
    show_output(f"EBcache for split {i} of {chrom} done! Gzipped cache saved to {cache_file}", color='success', time=True)
else:
    open(cache_file, 'a').close()
    show_output(f"Matrix file {matrix_file} was empty! Empty file {cache_file} is touched lest the pipeline break.", color='success', time=True)
