import pandas as pd
from HDR_utils import getHDR
from script_utils import show_output


def run_HDR(mut_file, tumor_bam, normal_bam, filter_pileup, out_file, MINSIM=.90, PAD=100, MINQ=25, HDRMINCOUNT=1):
    '''

    '''

    show_output(f'Starting HDR analysis of {mut_file}. [MIN_SIM={MINSIM}, PAD={PAD}]')
    # GET THE mutation file for masterHDR
    show_output(f'Importing {mut_file} for HDR detection', time=False)
    HDR_df = pd.read_csv(mut_file, sep='\t').loc[:, [
        'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']]

    # check for empty file
    if len(HDR_df.index) == 0:
        # return an empty csv with the right columns
        pd.DataFrame(columns=[
            'Chr',
            'Start',
            'End',
            'Ref',
            'Alt',
            'Gene',
            'TumorHDRcand',
            'TumorHDRcount',
            'TumorHDRinfo',
            'NormalHDRcand',
            'NormalHDRcount',
            'NormalHDRinfo'
        ]).to_csv(out_file, sep='\t', index=False)
        show_output(f"Mutation file {mut_file} is empty! Writing empty file {out_file}", color="warning")
        return

    # run HDR for both tumor and normal
    bam_dict = {'Tumor': tumor_bam, 'Normal': normal_bam}

    # ###### PILEUP ANALYSIS ##############################
    for T_or_N in bam_dict.keys():
        show_output(f'Analysing {T_or_N}')
        HDR_df = getHDR(bam_dict[T_or_N], HDR_df, _T_or_N=T_or_N, pileup_file=filter_pileup,
                        MINSIM=MINSIM, padding=PAD, min_HDR_count=HDRMINCOUNT, MINQ=MINQ)
        HDR_df = HDR_df.rename(columns={
            'HDR': f'{T_or_N}HDRcand',
            'HDRcount': f'{T_or_N}HDRcount',
            'HDRinfo': f'{T_or_N}HDRinfo'
        })   

    HDR_len = len(HDR_df.query('TumorHDRcount > 0').index)

    show_output(
        f"Found {HDR_len} possible HDR mutations. Writing to file {out_file}")
    HDR_df.to_csv(out_file, sep='\t', index=False)
s