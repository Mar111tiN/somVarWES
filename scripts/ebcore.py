import pandas as pd
import numpy as np
import math
from functools import partial
from scipy.optimize import fmin_l_bfgs_b as minimize_func
from scipy.stats import chi2
from scipy.special import gammaln


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


def get_target_count(row, var):
    '''
    converts the base-wise read coverage to a matrix
    '''

    s = pd.Series()

    s['depth_p'] = int(row['Depth-ACGT'])
    s['depth_n'] = int(row['Depth-acgt'])
    # treat Indels as combined event
    if var in ['I', 'D']:
        s['depth_p'] += int(row['Depth-INDEL'])
        s['depth_n'] += int(row['Depth-indel'])
        # sum up insert and del events
        s['mm_p'] = int(row['I']) + int(row['D'])
        s['mm_n'] = int(row['i']) + int(row['d'])
    else:
        s['mm_p'] = int(row[var])
        s['mm_n'] = int(row[var.lower()])
    return s


def matrix2EBscore(pen, row):

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


def AB2EBscore(row):
    '''
    get the EBscore from cached AB parameters
    no fitting is needed as parameters are precomputed and stored in row[5:9]
    '''

    # set the variant to the value fitting to the matrix
    if row['Ref'] == "-":
        ALT = "I"
    elif row['Alt'] == "-":
        ALT = "D"
    else:
        ALT = row['Alt'].upper()

    # we only get the snp count_df, using the mut_df 'Alt' as var and adjust for AB_df with column 9
    count_series = get_target_count(row, ALT)
    bb_params = {}

    # adjust variant D-->I and d-->i because I/i is used in ABcache
    if ALT == "D":
        ALT = ("I")
    # feed-in the AB params coming with the row
    bb_params['p'] = [row[f"{ALT}+a"], row[f"{ALT}+b"]]
    bb_params['n'] = [row[f"{ALT}-a"], row[f"{ALT}-b"]]

    p_values = bb_pvalues(bb_params, count_series)

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
        obs_list = [target_df + np.array([0, i])
                    for i in range(0, n_minus_k + 1)]
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
KS_matrix = np.array([[1, 0, 1, 1, 0, 1, 0, 0, 0], [
                     0, 1, -1, 0, 1, -1, 0, 0, 0]])
gamma_reduce = np.array([1, -1, -1, -1, 1, 1, 1, -1, -1])


def bb_loglikelihood(params, count_df, is_1d):
    [a, b] = params
    ab_matrix = np.array([1, 1, 1, a + b, a, b, a + b, a, b])
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
        result = penalty * \
            math.log(sum(params)) - bb_loglikelihood(params,
                                                     count_df, False)  # matrix is dim2
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


# ######################### EB-CACHE ####################

def get_cache_count_df(row, var):
    '''
    converts the base-wise read coverage to a matrix
    '''

    matrix = pd.DataFrame()
    matrix['depth_p'] = np.array(row['Depth-ACGT'].split('|')).astype(int)
    matrix['mm_p'] = np.array(row[var].split('|')).astype(int)
    matrix['depth_n'] = np.array(row['Depth-acgt'].split('|')).astype(int)
    matrix['mm_n'] = np.array(row[var.lower()].split('|')).astype(int)
    if var == 'I':
        matrix['depth_p'] += np.array(row['Depth-INDEL'].split('|')
                                      ).astype(int)
        matrix['depth_n'] += np.array(row['Depth-indel'].split('|')
                                      ).astype(int)
    return matrix


def matrix2AB(matrix_df, pen):
    '''
    creates the AB_df for a pileup
    '''

    # creates an empty-shell df from the matrix_df to receive the AB parameters
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
            count_df = get_cache_count_df(row, var)
            # get the AB parameters for

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
    var_columns = [f'{var}{strand}{param}' for var in list(
        'ACTGI') for strand in ['+', '-'] for param in ['a', 'b']]


# !!!!!!!!!!!!!!!!!! PEN!!!!!!!!!!!!!!!!!!!
    AB_df[var_columns] = matrix_df.apply(partial(get_AB, pen), axis=1)
    return AB_df
