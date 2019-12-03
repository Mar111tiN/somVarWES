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
######################## UTILS ############################

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

###########################################################
###################### EB MATH CORE #######################

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
 

###############################################################
###################### EBscore ROW HANDLERS ###################

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
    alt = ALT.lower()
    
    
    # we only get the snp count_df, using the mut_df 'Alt' as var and adjust for AB_df with column 9
    count_series = get_target_count(row, ALT)
    bb_params = {}
    
    # adjust variant D-->I and d-->i because I/i is used in ABcache
    if ALT == "D":
        ALT, alt = ("I", "i")
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


####################################################################
########################### SNAKE PARAMS ###########################

w = snakemake.wildcards
config = snakemake.config
threads = snakemake.threads
log = snakemake.log
output = snakemake.output
params = snakemake.params
pon_list = os.path.join(config['paths']['mystatic'], config['EBFilter']['pon_list'])
EBparams = snakemake.config['EBFilter']['params']
fit_pen = EBparams['fitting_penalty']
mut_file = snakemake.input.anno
tumor_bam = snakemake.input.tumor_bam
cache_file = snakemake.input.EBcache
matrix_cache = cache_file.replace('.cache', '.matrix')
chrom = w.chrom

# import the scripts
cleanpileup = params.cleanpileup
csv2bed = params.csv2bed
pon2cols = params.pon2cols
pile2count = params.pile2count
matrix2EBinput = params.matrix2EBinput
reduce_matrix = params.reducematrix
matrix2EBinput = params.matrix2EBinput


####################################################################
########################### RUN COMPUTATION ########################

############### LOAD DATA ###############################
show_output(f"Computing EBscore for chrom {chrom} of {tumor_bam} using EBcache {cache_file}", color='normal', time=True)

# get the mutation file for the chromosome
mut_df = pd.read_csv(mut_file, sep='\t', index_col=False).query('Chr == @chrom').iloc[:,:5]
mut_cols = list(mut_df.columns)

# set base for intermediate files
base_file = output[0].replace(".cachedEB","")


############### TARGET PILEUP --> MATRIX FILE ##################
mut_chr_file = f"{base_file}.csv"
mut_df.to_csv(mut_chr_file, sep='\t', index=False)
bed_file = f"{base_file}.bed"
# create the bed file for mpileup from the mutation file
shell(f"{csv2bed} < {mut_chr_file} > {bed_file}")

# do the pileup into the matrix file
tumor_matrix_file = f"{base_file}.matrix"
pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']} -l {bed_file} -r {chrom} {tumor_bam}"
pipe_cmd = f"{pileup_cmd} | cut -f 1,2,5 | {cleanpileup} | {pile2count} > {tumor_matrix_file}"
show_output(f"Piling up tumor bam {tumor_bam}", color='normal', time=True)
# do the pileup to matrix_file
show_command(pipe_cmd, multi=False)
shell(pipe_cmd)
# cleanup
shell(f"rm -f {bed_file} {mut_chr_file}")


show_output(f"Pileup matrix for chrom {chrom} of {tumor_bam} completed. Merging with cache file...", color='normal', time=True)

############### LOAD AND MERGE CACHE INTO MATRIX FILE #####

## check if sample is included in PoN ######
 def checkPon4sample(row, bam):
    '''
    per row of Pon_list checks whether sample name occurs in Pon_list
    if yes, returns the row number (zero-based)
    '''
    if os.path.basename(bam).split('_')[0] in row.iloc[0]: 
        print(row[0], row.name) 
        return row.name
    return 0

def get_sample_pos(pon_list, bam_file):
    '''
    returns the zero-based position of corresponding normal sample in Pon_list
    '''

    pon_df = pd.read(pon_list)
    sample_pos = pon_df.apply(checkPon4sample, axis=1, bam=tumor_bam).sum()
    return sample_pos

# if sample_inpon == 0, then sample is not in PoN
# else, pon matrix has to be acquired from cache and used in EBscore
sample_in_pon = get_sample_pos(pon_list, tumor_bam)


if sample_in_pon:
    ############################################ CACHE FROM MATRIX #####################################

    show_output(f"Corresponding normal sample for {tumor_bam} has been found in Panel of Normals. EBcache cannot be used!", color='warning')
    show_output(f"Falling back to cached matrix file reduced by corresponding normal..", color='normal')
    # get the cached matrix file reduced by the sample
    # EBcache cannot be used directly
    reduced_matrix_file = f"{base_file}.ponmatrix"
    reduce_matrix_cmd = f"gunzip < {matrix_file} | {reduce_matrix} {sample_in_pon} > {reduced_matrix_file}"

    ############# LOAD AND MERGE MATRIX FILES INTO MUTFILE

    # change mutation positions for deletions in mutation file
    mut_df.loc[mut_df['Alt'] == "-",'Start'] = mut_df['Start'] -1

    # load in the target matrix file as df
    tumor_matrix_df = pd.read_csv(tumor_matrix_file, sep='\t', index_col=False)
    # merge
    mut_matrix = mut_df.merge(tumor_matrix_df, on=['Chr', 'Start'], how='inner')
    # reset deletion positions
    mut_matrix.loc[mut_matrix['Alt'] == "-",'Start'] = mut_matrix['Start'] + 1

    # load in the pon_matrix file as df
    pon_matrix_df = pd.read_csv(reduced_matrix_file, sep='\t', index_col=False)
    # merge
    mut_matrix = mut_df.merge(pon_matrix_df, on=['Chr', 'Start'], how='inner')

    # write to file
    mutmatrix_file = f"{base_file}.mutmatrix"
    mut_matrix.to_csv(mutmatrix_file, sep='\t', index=False)


    ## CONTINUE LIKE UNCACHED EBscore
    # convert mutmatrix to direct EBinput
    EB_matrix_input_file = f"{base_file}.EB.matrix"
    shell(f"cat {mutmatrix_file} | {matrix2EBinput} > {EB_matrix_input_file}")

    # load in the EB.matrix file
    eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')

    # multithreaded computation
    def  computeEB2(df):
        show_output(f"Computing EBscore for {len(df.index)} lines", multi=True, time=True)
        df['EBscore'] = df.apply(partial(get_EB_score2, fit_pen), axis=1)
        show_output("Finished!", multi=True, time=True)
        return df

    eb_pool = Pool(threads)
    # minimal length of 1000 lines
    split_factor = min(math.ceil(len(eb_matrix.index) / 1000), threads)
    mut_split = np.array_split(eb_matrix, split_factor)
    EB_dfs = eb_pool.map(computeEB2, mut_split)
    
else:
    ############################################ CACHE FROM ABcache #####################################

    # get the mutation file for the chromosome
    cache_df = pd.read_csv(cache_file, compression='gzip', sep='\t')
    matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=False)
    pileAB_file = f"{base_file}.pileAB"
    pileAB_df = matrix_df.merge(cache_df, on=['Chr', 'Start'])
    # change coords for merge with start
    mut_df.loc[mut_df['Alt'] == "-",'Start'] = mut_df['Start'] -1
    pileAB_df = mut_df.merge(pileAB_df, on=['Chr', 'Start'])
    pileAB_df.loc[pileAB_df['Alt'] == "-",'Start'] = pileAB_df['Start'] + 1

    # save for debugging
    pileAB_df.to_csv(pileAB_file, sep='\t', index=False)
    show_output(f"Pileup matrix for for chrom {chrom} of {tumor_bam} merged with AB matrix.\n Written matrix file to {pileAB_file}.\n Going on with EB computation...", color='normal', time=True)


    ############### EBSCORE COMPUTATION  ########
    # multithreaded computation
    def  computeEB(df):
        show_output(f"Computing EBscore for {len(df.index)} lines", multi=True, time=True)
        df['EBscore'] = df.apply(AB2EBscore, axis=1)
        show_output("Finished!", multi=True, time=True)
        return df

    eb_pool = Pool(threads)
    # minimal length of 2000 lines
    split_factor = min(math.ceil(len(pileAB_df.index) / 2000), threads)
    pileAB_split = np.array_split(pileAB_df, split_factor)
    EB_dfs = eb_pool.map(computeEB, pileAB_split)



##### CONCAT AND WRITE TO FILE
EB_df = pd.concat(EB_dfs).sort_values(['Chr', 'Start'])

########## WRITE TO FILE ##############################################
# add EBscore to columns
mut_cols.append('EBscore')

# condense base info
base_cols = list("AaGgCcTtIiDd")
col_name = "|".join(base_cols)
# convert base coverage to str
for ch in base_cols:
    EB_df[ch] = EB_df[ch].map(str)
# condense base info into col "A|a|G|g|C|c|T|t|I|i|D|d"
EB_df[col_name] = EB_df[base_cols].apply(lambda row: "|".join(row), axis=1)
# add "A|a|G|g|C|c|T|t|I|i|D|d" to columns
mut_cols.append(col_name)
# rm unnecessary columns
EB_df = EB_df[mut_cols]

EB_df.to_csv(output[0], sep='\t', index=False)

shell(f"rm -f {matrix_file} {pileAB_file}") # 
show_output(f"Created EBscore for chrom {chrom} of {tumor_bam} using EBcache and written to {output[0]}", color='success', time=True)
