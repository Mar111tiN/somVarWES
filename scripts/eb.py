import os
import pandas as pd
from os import system as shell
from ebutils import get_pon_bases, compute_matrix2EB_multi
from script_utils import show_output, show_command

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
show_output(f"Computing EBscore for chrom {chrom} of {tumor_bam}", color='normal')

# get the sceleton mutation file
mut_df = pd.read_csv(mut_file, sep='\t', index_col=False).query('Chr == @chrom').iloc[:, :5]
mut_cols = list(mut_df.columns)
# set base_name for intermediate files
base_file = output[0].replace(".EB", "")

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

show_output(f"Piling up {chrom} of {tumor_bam} with Pon List.", color='normal')
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
    show_output(f"Pileup matrix for chrom {chrom} of {tumor_bam} completed.", color='normal')
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
    EB_df = compute_matrix2EB_multi(eb_matrix, fit_pen, threads)

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
        print('full_output')
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
