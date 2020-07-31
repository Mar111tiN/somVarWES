import os
import pandas as pd
from os import system as shell
from ebutils import get_pon_bases, get_sample_pos, compute_matrix2EB_multi, compute_AB2EB_multi
from script_utils import show_output, show_command, run_cmd

def run_eb(table, tumor_bam, output, pon_list, chrom, log, threads, EBparams, full_output,
           cleanpileup,
           csv2bed,
           pon2cols,
           pile2count,
           matrix2EBinput,
           makeponlist
    ):
    '''
    master function to start eb_computation
    '''
    # ############## LOAD DATA ###############################
    show_output(f"Computing EBscore for chrom", color='normal')

    # get the sceleton mutation file for that chromosome
    mut_df = pd.read_csv(table, sep='\t', index_col=False, header=None, names=['Chr', 'Start', 'End', 'Ref', 'Alt', 'somatic_status', 'TR1', 'TR1+', 'TR2', 'TR2+', 'NR1', 'NR1+', 'NR2', 'NR2+', 'somaticP', 'variantP']).query('Chr == @chrom').iloc[:, :5]
    mut_cols = list(mut_df.columns)
    # set base_name for intermediate files
    base_file = output[0].replace(".EB", "")

    # ############## PILEUP --> MATRIX FILE ##################

    # bed file can contain all chromosomes because chrom restriction comes with the -r parameter
    bed_file = f"{base_file}.bed"
    # create the bed file for mpileup
    shell(f"{csv2bed} < {table} > {bed_file}")

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
    pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']}"
    pileup_cmd += f" -l {bed_file} -r {chrom} -b {sample_list}"
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
        mut_matrix.loc[mut_matrix['Alt'] == "-", 'Start'] = mut_matrix['Start'] + 1

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
        EB_df = compute_matrix2EB_multi(eb_matrix, EBparams['fitting_penalty'], threads)

        # add EBscore to columns
        mut_cols.append('EBscore')

        # get the pon_matrix containing the Pon coverages in Alt and Ref
        pon_matrix = get_pon_bases(eb_matrix)
        # transfer PoN-Ref and PoN-Alt to EB_df
        EB_df[['PoN-Ref', 'PoN-Alt']] = pon_matrix[['PoN-Ref', 'PoN-Alt']]
        mut_cols += ['PoN-Ref', 'PoN-Alt']

        # ###### add the full output ##########
        if full_output:
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
        shell(f"rm {matrix_file} {EB_matrix_input_file}")  # {mutmatrix_file}
        show_output(f"Created EBscore for chrom {chrom} of {tumor_bam} and written to {output[0]}", color='success')



def run_eb_from_cache(table, tumor_bam, output, pon_list, chrom, log, threads, EBparams, full_output,
           cleanpileup,
           csv2bed,
           pile2count,
           matrix2EBinput,
           reorder_matrix
    ):

    
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
        pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']}"
        pileup_cmd += f" -l {bed_file} -r {chrom} {tumor_bam}"
        pipe_cmd = f"{pileup_cmd} | cut -f 1,2,5 | {cleanpileup} | {pile2count} > {matrix_file}"
        show_output(f"Piling up tumor bam {tumor_bam}", color='normal')
        # do the pileup to matrix_file
        show_command(pipe_cmd, multi=False)
        shell(pipe_cmd)
        # cleanup
        shell(f"rm -f {bed_file}")
        show_output(
            f"Pileup matrix for chrom {chrom} of {tumor_bam} completed. Merging with cache file...",
            color='normal'
        )
        return matrix_file

    # ############## LOAD DATA ###############################
    show_output(f"Computing EBscore for chrom {chrom} of {tumor_bam} using EBcache {AB_cache_file}", color='normal')

    # get the mutation file for the chromosome
    mut_df = pd.read_csv(mut_file, sep='\t', index_col=False, header=None, names=['Chr', 'Start', 'End', 'Ref', 'Alt', 'somatic_status', 'TR1', 'TR1+', 'TR2', 'TR2+', 'NR1', 'NR1+', 'NR2', 'NR2+', 'somaticP', 'variantP']).query('Chr == @chrom').iloc[:, :5]
    mut_cols = list(mut_df.columns)

    # check for empty df
    if mut_df.empty:
        EB_df = pd.DataFrame(columns=mut_cols)
        EB_df.to_csv(output[0], sep='\t', index=False)
        show_output(f"No mutations for {chrom} in mutation list! Writing empty file to {output[0]}", color='warning')
    else:
        # set base_name for intermediate files
        base_file = output[0].replace(".cachedEB", "")

        # ############## LOAD PILEUP MATRIX CACHE AND MERGE INTO MUT_DF #####
        # change mutation positions for deletions in mutation file
        mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] - 1
        show_output(f"Loading compressed matrix cache file {matrix_cache_file}", color='normal')
        # load in the target matrix file as df
        cache_matrix_df = pd.read_csv(matrix_cache_file, sep='\t', index_col=False, compression='gzip')
        # merge
        mut_matrix = mut_df.merge(cache_matrix_df, on=['Chr', 'Start'], how='inner')
        # reset deletion positions
        mut_matrix.loc[mut_matrix['Alt'] == "-", 'Start'] = mut_matrix['Start'] + 1
        show_output(f"Loaded and merged into mutation list", color='normal')

        # ############### CHECK IF SAMPLE IN PON ####################
        # if sample_inpon == 0, then sample is not in PoN
        # else, pon matrix has to be acquired from cache and used in EBscore
        sample_in_pon = get_sample_pos(pon_list, tumor_bam)

        # ########################################### CACHE FROM MATRIX #####################################
        if sample_in_pon:
            show_output(
                f"Corresponding normal sample for {tumor_bam} has been found in PoNs! EBcache cannot be used!",
                color='warning'
            )
            show_output(f"Falling back to cached matrix file..", color='normal')
            # EBcache cannot be used directly

            # ######### REMOVE SAMPLE BASES FROM MATRIX FILE

            # get the cached matrix file and reorder sample bases to first position to create valid mutmatrix
            # reorder_matrix takes position of sample in pon_list as argument
            # if position of tumor bam in pon == 1, everything is already fine
            mutmatrix_file = f"{base_file}.mutmatrix"
            if sample_in_pon > 1:
                prematrix_file = f"{base_file}.prematrix"
                mut_matrix.to_csv(prematrix_file, sep='\t', index=False)

                # row is 0-based --> sample_in_pon + 1
                reduce_matrix_cmd = f"cat {prematrix_file} | {reorder_matrix} {sample_in_pon - 1} > {mutmatrix_file}"
                show_command(reduce_matrix_cmd, multi=False)
                shell(reduce_matrix_cmd)
                # cleanup
                shell(f"rm {prematrix_file}")
            else:
                # tumor sample already in the right position
                mut_matrix.to_csv(mutmatrix_file, sep='\t', index=False)
            show_output(f"Retrieving target data from cached matrix", color='normal')

            # # CONTINUE LIKE UNCACHED EBscore
            # convert mutmatrix to direct EBinput
            EB_matrix_input_file = f"{base_file}.EB.matrix"
            EBinput_cmd = f"cat {mutmatrix_file} | {matrix2EBinput} > {EB_matrix_input_file}"
            show_command(EBinput_cmd, multi=False)
            shell(EBinput_cmd)
            # load in the EB.matrix file
            eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')
            print('Start computation file')
            # multithreaded computation
            # passing attempts to threads
            EB_df = compute_matrix2EB_multi(eb_matrix, EBparams['fitting_penalty'], threads)
            print('Computation finished')
            # get the pon_matrix containing the Pon coverages in Alt and Ref
            pon_matrix = get_pon_bases(eb_matrix)
        # ########################################### CACHE FROM ABcache ###########################
        else:
            # ############## TARGET PILEUP --> MATRIX FILE ##################
            tumor_matrix_file = target_pileup_from_mut(mut_file, base_file, tumor_bam, chrom)
            # check if matrix_file has input
            # if not os.path.getsize(tumor_matrix_file):
            #     # create empty file
            #     EB_df = mut_df
            #     EB_df['EBscore'] = 0
            #     has_pileup = False

            # else:  # has input
            # has_pileup = True
            # reloading the target pileup into pileup_df
            # use dtype to ensure str encoding of chromosome columns
            pileup_df = pd.read_csv(tumor_matrix_file, sep='\t', dtype={'Chr': str, 'Start': int}, index_col=False)

            show_output(f"Loading compressed AB cache file {AB_cache_file}", color='normal')
            cache_df = pd.read_csv(AB_cache_file, compression='gzip', sep='\t')

            pileAB_df = pileup_df.merge(cache_df, on=['Chr', 'Start'])
            # change coords for merge with start and merge into mut_df for Ref
            mut_df.loc[mut_df['Alt'] == "-", 'Start'] = mut_df['Start'] - 1
            pileAB_df = mut_df.merge(pileAB_df, on=['Chr', 'Start'])
            pileAB_df.loc[pileAB_df['Alt'] == "-", 'Start'] = pileAB_df['Start'] + 1

            # save for debugging
            # pileAB_file = f"{base_file}.pileAB"
            # pileAB_df.to_csv(pileAB_file, sep='\t', index=False)
            show_output(
                f"Pileup matrix for for chrom {chrom} of {tumor_bam} merged with AB matrix." +
                " Going on with EB computation...",
                color='normal'
            )

            # ############## EBSCORE COMPUTATION  ########
            # multithreaded computation
            EB_df = compute_AB2EB_multi(pileAB_df, threads)

            # convert matrix file to EB_input for getting PoN-Ref and Pon-Alt
            mutmatrix_file = f"{base_file}.mutmatrix"
            mut_matrix.to_csv(mutmatrix_file, sep='\t', index=False)
            # do the conversion
            EB_matrix_input_file = f"{base_file}.EB.matrix"
            convert_cmd = (f"cat {mutmatrix_file} | {matrix2EBinput} > {EB_matrix_input_file}")
            show_command(convert_cmd)
            shell(convert_cmd)

            # load in the EB.matrix file
            eb_matrix = pd.read_csv(EB_matrix_input_file, sep='\t')

            # get the pon_matrix containing the Pon coverages in Alt and Ref
            # tumor sample is not in PoN --> no removal neccessary
            pon_matrix = get_pon_bases(eb_matrix, remove_sample=False)

            # cleanup
            shell(f"rm -f {tumor_matrix_file}")

        # add EBscore to columns
        mut_cols.append('EBscore')

        # transfer PoN-Ref and PoN-Alt from pon_matrix to EB_df
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
        EB_df.to_csv(output[0], sep='\t', index=False)

        # cleanup
        shell(f"rm -f {EB_matrix_input_file}")
        show_output(
            f"Created EBscore for chrom {chrom} of {tumor_bam} using EBcache and written to {output[0]}",
            color='success'
        )