import os
import pandas as pd
from script_utils import show_output, sort_df


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    i = s.input
    o = s.output
    p = s.params
    c = s.config
    cc = c["HDR"]

    # populate an input dict for better looping
    input_list = {'Tumor': i.tumor, 'Normal': i.normal}
    base_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt']
    cols = base_cols.copy()
    dfs = {}
    for T_or_N in ['Tumor', 'Normal']:
        HDR_dfs = []
        for HDR_file in input_list[T_or_N]:
            if os.path.isfile(HDR_file):
                try:
                    HDR_df = pd.read_csv(HDR_file, sep='\t', index_col=False)
                    if HDR_df.empty:
                        continue
                    HDR_dfs.append(HDR_df)
                except:
                    show_output(f"{HDR_file} could not be loaded")
        # concat
        if len(HDR_dfs):
            HDR_df = pd.concat(HDR_dfs).sort_values(['Chr', 'Start'])
        else:
            HDR_df = pd.DataFrame(columns=base_cols + ['HDRcand', 'HDRcount', 'HDRinfo'])
        # change columns
        new_cols = {
            'HDRcand': f'{T_or_N}HDRcand',
            'HDRcount': f'{T_or_N}HDRcount',
            'HDRinfo': f'{T_or_N}HDRinfo'
        }

        HDR_df = HDR_df.rename(columns=new_cols)
        # add new_cols to cols
        cols += new_cols.values()
        # load into dfs dict
        dfs[T_or_N] = HDR_df
    
    HDR_merge = pd.merge(dfs['Tumor'], dfs['Normal'], on=base_cols)
    # sort columns
    HDR_merge = sort_df(HDR_merge[cols])

    # ###### PILEUP ANALYSIS ##############################

    HDR_len = len(HDR_merge.query('TumorHDRcount > 0').index)
    HDR_merge.to_csv(str(o), sep='\t', index=False)
    show_output(f"Found {HDR_len} possible HDR mutations. Written HDR file to {str(o)}", color='success')


if __name__ == "__main__":
    main(snakemake)
