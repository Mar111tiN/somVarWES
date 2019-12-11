import os
from os import system as shell
import pandas as pd
from ebutils import show_output, show_command, get_pon_bases, get_sample_pos, compute_matrix2EB_multithreaded, computeEBfromAB_multithreaded


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
AB_cache = snakemake.input.EBcache
matrix_cache = AB_cache.replace('.cache', '.matrix')
chrom = w.chrom

# import the scripts
cleanpileup = params.cleanpileup
csv2bed = params.csv2bed
pon2cols = params.pon2cols
pile2count = params.pile2count
matrix2EBinput = params.matrix2EBinput
reduce_matrix = params.reducematrix
matrix2EBinput = params.matrix2EBinput
reorder_matrix = params.reorder_matrix


def target_pileup_from_mut(mut_file, base_file, bam, chrom):
    '''
    piles up the mutation list in the tumor bam
    '''

    # bed file can contain all chromosomes because chrom restriction comes with the -r parameter
    bed_file = f"{base_file}.bed"
    # create the bed file for mpileup
    shell(f"{csv2bed} < {mut_file} > {bed_file}")

    # # if I want to restrict chromosome in file:
    # mut_chr_file = f"{base_file}.csv"
    # mut_df.to_csv(mut_chr_file, sep='\t', index=False)
    # # create the bed file for mpileup from the mutation file
    # shell(f"{csv2bed} < {mut_chr_file} > {bed_file}")

    # do the pileup into the matrix file
    matrix_file = f"{base_file}.matrix"
    pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']} -l {bed_file} -r {chrom} {tumor_bam}"
    pipe_cmd = f"{pileup_cmd} | cut -f 1,2,5 | {cleanpileup} | {pile2count} > {matrix_file}"
    show_output(f"Piling up tumor bam {tumor_bam}", color='normal', time=True)
    # do the pileup to matrix_file
    show_command(pipe_cmd, multi=False)
    shell(pipe_cmd)
    # cleanup
    shell(f"rm -f {bed_file}")
    show_output(f"Pileup matrix for chrom {chrom} of {tumor_bam} completed. Merging with cache file...", color='normal', time=True)
    return matrix_file


# ############## LOAD DATA ###############################
show_output(f"Computing EBscore for chrom {chrom} of {tumor_bam} using EBcache {AB_cache}", color='normal', time=True)

# get the mutation file for the chromosome
mut_df = pd.read_csv(mut_file, sep='\t', index_col=False).query('Chr == @chrom').iloc[:, :5]
mut_cols = list(mut_df.columns)
# set base_name for intermediate files
base_file = output[0].replace(".cachedEB","")

# ############## LOAD PILEUP MATRIX CACHE AND MERGE INTO MUT_DF #####
# change mutation positions for deletions in mutation file
mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] -1
# load in the target matrix file as df
tumor_matrix_df = pd.read_csv(matrix_cache, sep='\t', index_col=False)
# merge
mut_matrix = mut_df.merge(tumor_matrix_df, on=['Chr', 'Start'], how='inner')
# reset deletion positions
mut_matrix.loc[mut_matrix['Alt'] == "-", 'Start'] = mut_matrix['Start'] + 1


# ############### CHECK IF SAMPLE IN PON ####################
# if sample_inpon == 0, then sample is not in PoN
# else, pon matrix has to be acquired from cache and used in EBscore
sample_in_pon = get_sample_pos(pon_list, tumor_bam)


# ########################################### CACHE FROM MATRIX #####################################
if sample_in_pon:
    show_output(f"Corresponding normal sample for {tumor_bam} has been found in Panel of Normals. EBcache cannot be used!", color='warning')
    show_output(f"Falling back to cached matrix file..", color='normal')
    # EBcache cannot be used directly

    # ######### REMOVE SAMPLE BASES FROM MATRIX FILE
    # get the cached matrix file and reorder sample bases to first position to create valid mutmatrix
    # reorder_matrix takes position of sample in pon_list as argument
    # write to file
    prematrix_file = f"{base_file}.prematrix"
    mutmatrix_file = f"{base_file}.mutmatrix"
    mut_matrix.to_csv(prematrix_file, sep='/t', index=False)
    reduce_matrix_cmd = f"cat {prematrix_file} | {reorder_matrix} {sample_in_pon} > {mutmatrix_file}"
    # reload reordered mutmatrix into mut_matrix
    mut_matrix = pd.read_csv(mutmatrix_file, sep='\t', index_col=False)
    # cleanup
    shell(f"rm {prematrix_file}")

    # # CONTINUE LIKE UNCACHED EBscore
    # convert mutmatrix to direct EBinput
    EB_matrix_input_file = f"{base_file}.EB.matrix"
    shell(f"cat {mutmatrix_file} | {matrix2EBinput} > {EB_matrix_input_file}")

    # load in the EB.matrix file
    eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')

    # multithreaded computation
    EB_df = compute_matrix2EB_multithreaded(eb_matrix, fit_pen, threads)
    # add EBscore to columns
    mut_cols.append('EBscore')

    # get the pon_matrix containing the Pon coverages in Alt and Ref
    pon_matrix = get_pon_bases(eb_matrix)
    # transfer PoN-Ref and PoN-Alt to EB_df
    EB_df[['PoN-Ref', 'PoN-Alt']] = pon_matrix[['PoN-Ref', 'PoN-Alt']]
    mut_cols += ['PoN-Ref', 'PoN-Alt']

# ########################################### CACHE FROM ABcache ###########################
else:
    # ############## TARGET PILEUP --> MATRIX FILE ##################
    tumor_matrix_file = target_pileup_from_mut(mut_file, base_file, tumor_bam, chrom)



    # check if matrix_file has input
    if not os.path.getsize(tumor_matrix_file):
        # create empty file
        open(output[0], 'a').close()
        show_output(f"Pileup for {chrom} of {tumor_bam} was empty! Created empty file {output[0]}", color='warning')
    else:  # has input
        # get the mutation file for the chromosome
        cache_df = pd.read_csv(AB_cache, compression='gzip', sep='\t')
        matrix_df = pd.read_csv(matrix_cache, sep='\t', index_col=False)
        pileAB_file = f"{base_file}.pileAB"
        pileAB_df = matrix_df.merge(cache_df, on=['Chr', 'Start'])
        # change coords for merge with start
        mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] -1
        pileAB_df = mut_df.merge(pileAB_df, on=['Chr', 'Start'])
        pileAB_df.loc[pileAB_df['Alt'] == "-", 'Start'] = pileAB_df['Start'] + 1

        # save for debugging
        pileAB_df.to_csv(pileAB_file, sep='\t', index=False)
        show_output(f"Pileup matrix for for chrom {chrom} of {tumor_bam} merged with AB matrix.\n Written matrix file to {pileAB_file}.\n Going on with EB computation...", color='normal', time=True)

        # ############## EBSCORE COMPUTATION  ########
        # multithreaded computation
        EBdf = computeEBfromAB_multithreaded(pileAB_df, fit_pen, threads)

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

    # ######### WRITE TO FILE ##############################################
    EB_df.to_csv(output[0], sep='\t', index=False)

    # cleanup
    shell(f"rm -f {matrix_file} {pileAB_file}") # 
    show_output(f"Created EBscore for chrom {chrom} of {tumor_bam} using EBcache and written to {output[0]}", color='success', time=True)
