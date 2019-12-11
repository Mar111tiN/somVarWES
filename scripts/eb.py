import os
from os import system as shell
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
# ##################### FUNCTIONS ##########################

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


# ############## EBscore using matrix2EBinput.mawk #######################
def get_count_df(row):
    '''
    converts the base-wise read coverage to a matrix
    '''

    matrix = pd.DataFrame()
    matrix['depth_p'] = np.array(row['depthP'].split('|')).astype(int)
    matrix['mm_p'] = np.array(row['misP'].split('|')).astype(int)
    matrix['depth_n'] = np.array(row['depthN'].split('|')).astype(int)
    matrix['mm_n'] = np.array(row['misN'].split('|')).astype(int)
    return matrix


def get_EB_score2(pen, row):

    count_df = get_count_df(row)
    # ########### FITTING ####################################
    # get the respective control matrices (as dataframe) for positive and negative strands
    # estimate the beta-binomial parameters for positive and negative strands

    # <<<<<<######### DEBUG ###############
    # print(row['Chr'], row['pos'], count_df)
    # <<<<<<###############################

    bb_params = fit_bb(count_df[1:], pen)
    # evaluate the p-values of target mismatch numbers for positive and negative strands
    p_values = bb_pvalues(bb_params, count_df.iloc[0])

    # ########### FISHER COMBINATION #########################
    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = fisher_combination(p_values)
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = -round(math.log10(EB_pvalue), 3)
    return EB_score
############################################################################


def fisher_combination(p_values):

    if 0 in p_values.values():
        return 0
    else:
        return 1 - chi2.cdf(sum([-2 * math.log(x) for x in p_values.values()]), 2 * len(p_values.values()))


def bb_pvalues(params, target_df):
    '''
    accumulate p_value of target observation falling in fitted bb_distribution (not a variant)
    p_values are computed per strand (pvalue_p and pvalue_n)
    p_value: exponential sum of loglikelihooks of successes greater or equal than observed
    [n, k] --> sum of density (exp of loglikelihood) [n, k] to [n, n]
    '''

    def bb_pvalue(params, target_df):
        n_minus_k = target_df[0] - target_df[1]
        # get the list of observations [n, k] to [n, n]
        obs_list = [target_df + np.array([0,i]) for i in range(0, n_minus_k + 1)]
        # get the list of loglikelihoods per observation
        ll_list = [bb_loglikelihood(params, obs, True) for obs in obs_list]

        #######################################################
        # print(f'ab: {params}\n observations: {obs_list} ll {ll_list}\n')
        #######################################################

        # get the sum of exponentials of loglikelihoods (densities) per observation

        p_value = sum([math.exp(ll) for ll in ll_list])

        return p_value

    target_p = target_df.loc[['depth_p', 'mm_p']]
    target_n = target_df.loc[['depth_n', 'mm_n']]
    p_values = {}
    p_values['p'] = bb_pvalue(params['p'], target_p)
    p_values['n'] = bb_pvalue(params['n'], target_n)
    return p_values


# the matrices for beta-binomial calculation
KS_matrix = np.array([[1, 0, 1, 1, 0, 1, 0, 0, 0], [0, 1, -1, 0, 1, -1, 0, 0, 0]])
gamma_reduce = np.array([1, -1, -1, -1, 1, 1, 1, -1, -1])


def bb_loglikelihood(params, count_df, is_1d):
    [a, b] = params
    ab_matrix = np.array([1, 1, 1, a+b, a, b, a+b, a, b])
    # convert df into matrix for np.array operations that change dims
    count_matrix = count_df.values
    # perform matrix multiplication to get inputs to log-gamma
    input_matrix = np.matmul(count_matrix, KS_matrix) + ab_matrix
    # get corresponding log-gamma values and reduce over pon-values
    if is_1d:  # check whether gammatrix is 2-dim - otherwise sum aggregation over axis 0 is faulty
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
        result = penalty * math.log(sum(params)) - bb_loglikelihood(params, count_df, False)  # matrix is dim2
        return result

    # get the respective control matrices (as dataframe) for positive and negative strands
    count_p = count_df.loc[:, ['depth_p', 'mm_p']]
    count_n = count_df.loc[:, ['depth_n', 'mm_n']]
    # minimize loglikelihood using L-BFGS-B algorithm
    ab_p = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args=(count_p, pen), approx_grad=True,
                           bounds=[(0.1, 10000000), (1, 10000000)]
                          )[0]
    ab_p = [round(param, 5) for param in ab_p]
    ab_n = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args=(count_n, pen), approx_grad=True,
                           bounds=[(0.1, 10000000), (1, 10000000)]
                          )[0]
    ab_n = [round(param, 5) for param in ab_n]
    return {'p': ab_p, 'n': ab_n}


# ######### ADD BASE INFO ##############################################

def get_pon_bases(matrix_df):
    '''
    returns from eb-matrix file the concatenated Pon coverage for pos and neg strand
    this is important output for mutation QC
    imput cols:
        depthP
        depthN
        misP
        misN
    '''

    # remove sample depths from the columns
    for col in ['depthP', 'misP', 'depthN', 'misN']:
        matrix_df[col] = matrix_df[col].str.replace(r"^[0-9]+\|","")

    # concate the respective columns
    matrix_df['PoN-Ref'] = matrix_df['depthP'].str.cat(matrix_df['depthN'], sep="-")
    matrix_df['PoN-Alt'] = matrix_df['misP'].str.cat(matrix_df['misN'], sep="-")
    return matrix_df
#######################################################################
#######################################################################


w = snakemake.wildcards
config = snakemake.config
threads = snakemake.threads
log = snakemake.log
output = snakemake.output
params = snakemake.params
pon_list = os.path.join(config['paths']['mystatic'], config['EBFilter']['pon_list'])
EBparams = snakemake.config['EBFilter']['params']
fit_pen = EBparams['fitting_penalty']
mut_file = snakemake.input.table
tumor_bam = snakemake.input.tumor_bam
chrom = w.chrom

# import the scripts
cleanpileup = params.cleanpileup
csv2bed = params.csv2bed
pon2cols = params.pon2cols
pile2count = params.pile2count
matrix2EBinput = params.matrix2EBinput
makeponlist = params.makeponlist

# ############## LOAD DATA ###############################
show_output(f"Computing EBscore for chrom {chrom} of {tumor_bam}", color='normal', time=True)

# get the sceleton mutation file
mut_df = pd.read_csv(mut_file, sep='\t', index_col=False).query('Chr == @chrom').iloc[:,:5]
mut_cols = list(mut_df.columns)
# set base_name for intermediate files
base_file = output[0].replace(".EB","")

# ############## PILEUP --> MATRIX FILE ##################

# bed file can contain all chromosomes because chrom restriction comes with the -r parameter
bed_file = f"{base_file}.bed"
# create the bed file for mpileup
shell(f"{csv2bed} < {mut_file} > {bed_file}")

# # if I want to restrict chromosome in file:
# mut_chr_file = f"{base_file}.csv"
# mut_df.to_csv(mut_chr_file, sep='\t', index=False)
# # create the bed file for mpileup from the mutation file
# shell(f"{csv2bed} < {mut_chr_file} > {bed_file}")


# create the pon_list containing the tumor-bam as first file
sample_list = f"{base_file}.pon"
# makeponlist removes the sample itself from list if it is part of PoN
shell(f"{makeponlist} {tumor_bam} {pon_list} {sample_list}")

show_output(f"Piling up {chrom} of {tumor_bam} with Pon List.", color='normal', time=True)
shell(f"cat {sample_list}")
# do the pileup into the matrix file
matrix_file = f"{base_file}.matrix"
pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']} -l {bed_file} -r {chrom} -b {sample_list}"
# cut -f $({pon2cols}< {sample_list}) creates a cut command only including the desired

pipe_cmd = f"{pileup_cmd} | cut -f $({pon2cols} < {sample_list}) | {cleanpileup} | {pile2count} > {matrix_file}"
# do the pileup to matrix_file
show_command(pipe_cmd, multi=False)
shell(pipe_cmd)
# cleanup
shell(f"rm {bed_file} {sample_list}")

# check if matrix_file has input
if not os.path.getsize(matrix_file):
    # create empty file
    open(output[0], 'a').close()
    show_output(f"Pileup for {chrom} of {tumor_bam} was empty! Created empty file {output[0]}", color='warning')
else:
    show_output(f"Pileup matrix for chrom {chrom} of {tumor_bam} completed.", color='normal', time=True)
    # ################ MERGE INTO MUTFILE ######################
    # change mutation positions for deletions in mutation file
    mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] - 1
    # read matrix file into df
    matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=False)
    # merge
    mut_matrix = mut_df.merge(matrix_df, on=['Chr', 'Start'], how='inner')
    # reset deletion positions
    mut_matrix.loc[mut_matrix['Alt'] == "-",'Start'] = mut_matrix['Start'] + 1

    # ####### if using matrix2EBinput.mawk #######################
    # write to file
    mutmatrix_file = f"{base_file}.mutmatrix"
    mut_matrix.to_csv(mutmatrix_file, sep='\t', index=False)

    # convert mutmatrix to direct EBinput
    EB_matrix_input_file = f"{base_file}.EB.matrix"
    shell(f"cat {mutmatrix_file} | {matrix2EBinput} > {EB_matrix_input_file}")
    # load in the EB.matrix file
    eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')

    # multithreaded computation
    def computeEB2(df):
        show_output(f"Computing EBscore for {len(df.index)} lines", multi=True, time=True)
        df['EBscore'] = df.apply(partial(get_EB_score2, fit_pen), axis=1)
        show_output("Finished!", multi=True, time=True)
        return df

    eb_pool = Pool(threads)
    # minimal length of 1000 lines
    split_factor = min(math.ceil(len(eb_matrix.index) / 1000), threads)
    mut_split = np.array_split(eb_matrix, split_factor)
    EB_dfs = eb_pool.map(computeEB2, mut_split)

    # EB_df contains EB_score and matrix for Ref and Alt
    EB_df = pd.concat(EB_dfs).sort_values(['Chr', 'Start'])

    # add EBscore to columns
    mut_cols.append('EBscore')

    # get the pon_matrix containing the Pon coverages in Alt and Ref
    pon_matrix = get_pon_bases(eb_matrix)
    # transfer PoN-Ref and PoN-Alt to EB_df
    EB_df[['PoN-Ref', 'PoN-Alt']] = pon_matrix[['PoN-Ref', 'PoN-Alt']]
    mut_cols += ['PoN-Ref', 'PoN-Alt']

    # ###### add the full output ##########
    if config['EBFilter']['full_pon_output']:
        # condense base info
        base_cols = list("AaGgCcTtIiDd")
        col_name = "|".join(base_cols)
        # convert base coverage to str
        for ch in base_cols:
            # take the letter info from the mut_matrix which is not yet condensated
            # str.replace removes the tumor bases
            EB_df[ch] = mut_matrix[ch].map(str).str.replace(r'^[0-9]+\|', "")
        # condense base info into col "A|a|G|g|C|c|T|t|I|i|D|d"
        EB_df[col_name] = EB_df[base_cols].apply(lambda row: "-".join(row), axis=1)
        # add "A|a|G|g|C|c|T|t|I|i|D|d" to columns
        mut_cols.append(col_name)
    # rm unnecessary columns
    EB_df = EB_df[mut_cols]

    # ######### WRITE TO FILE ##############################################

    EB_file = output[0]
    EB_df.to_csv(EB_file, sep='\t', index=False)

    # cleanup
    shell(f"rm {matrix_file} {EB_matrix_input_file}") # {mutmatrix_file} 
    show_output(f"Created EBscore for chrom {chrom} of {tumor_bam} and written to {output[0]}", color='success')
