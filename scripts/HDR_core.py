import pandas as pd
from HDR_utils import get_HDR_multi, get_filter_hdr_multi, get_hotspot_df
from script_utils import show_output


def HDR_master(mut_file, bam_file, chrom, pileup_file, threads, HDR_config, out_file=''):
    '''
    does all the work and gets directly called by HDR.py
    '''
    COLS = [
        'Chr',
        'Start',
        'End',
        'Ref',
        'Alt',
        'Gene',
        'HDRcand',
        'HDRcount',
        'HDRinfo',
    ]

    show_output(
        f"Starting HDR analysis of {bam_file} on chrom {chrom}. [MIN_SIM={HDR_config['MINSIM']}, PAD={HDR_config['PAD']}]")
    # GET THE mutation file for masterHDR
    show_output(f'Importing {mut_file} for HDR detection', time=False)
    mut_df = pd.read_csv(mut_file, sep='\t').loc[:, [
        'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']].query('Chr == @chrom')

    # check for empty file
    if len(mut_df.index) == 0:
        # return an empty df
        message = f"Mutation file {mut_file} for chrom {chrom} is empty!"
        if out_file:
            message += f" Writing empty file {out_file}"
            HDR_df = pd.DataFrame(columns=COLS)
            HDR_df.to_csv(out_file, sep='\t', index=False)
        show_output(message, color="warning")
        return HDR_df

    # get the right pileup_df (different cols for Tumor and Normal)
    hotspot_df = get_hotspot_df(pileup_file)
    show_output(
        f"Detected {len(hotspot_df.index)} putative HDR lanes in {bam_file}.")
    filter_HDR = get_filter_hdr_multi(mut_df, hotspot_df, threads, HDR_config)
    show_output(
        f"Found a total of {len(filter_HDR.index)} HDR-rich mutations", multi=False)

    HDR_df = get_HDR_multi(filter_HDR, bam_file=bam_file, hotspot_df=hotspot_df, threads=threads,
                           chrom=chrom, HDR_config=HDR_config)

    HDR_df = mut_df.merge(
        HDR_df, on=COLS[:6], how="outer").sort_values(COLS[:2])
    for col in COLS[6:8]:
        HDR_df[col] = HDR_df[col].fillna(0).astype(int)
    HDR_df['HDRinfo'] = HDR_df['HDRinfo'].fillna('no HDR in vincinity')
    if out_file:
        show_output(
            f"Finished HDR analysis of {bam_file} on chrom {chr}.", color='success')
        HDR_df.to_csv(out_file, sep='\t', index=False)
    return HDR_df
