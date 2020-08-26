import os
import pandas as pd
from script_utils import show_output


def sort_df(df, cols={'Chr': True, 'Start': True}):
    '''
    helper for sorting dfs for chromosomes using Chr, Start + cols in cols
    '''
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ['chrX', 'chrY']

    df['Chr'] = pd.Categorical(df['Chr'], chrom_list)
    return df.sort_values(list(cols.keys()), ascending=list(cols.values()))


def filter1(mut_file, basic_output, filter1_output,
            filter_file, filter_sheet, filter_name, keep_syn=False):

    show_output(f"Running filter1 \"{filter_name}\"")
    show_output(f"Loading mutation file {mut_file}.")
    anno_df = pd.read_csv(mut_file, sep='\t')

    show_output(f"Loading filter file {filter_file}")
    if "xls" in os.path.splitext(filter_file)[1]:
        filter_settings = pd.read_excel(
            filter_file, sheet_name=filter_sheet, index_col=0)[:4]
    else:
        filter_settings = pd.read_csv(filter_file, sep='\t', index_col=0)

    # use these population columns for checking PopFreq
    # could be refactored into params
    pop_cols = ['gnomAD_exome_ALL', 'esp6500siv2_all', 'dbSNP153_AltFreq']

    show_output(f"    keep_syn= {keep_syn}")
    show_output(f'Started editing and basic filtering for {mut_file}.')

    #  ############## BASIC FILTER ####################################
    def filter_basic(df, keep_syn=False):
        '''
        basic cutoff based on gene function
        '''

        exon_func = df['ExonicFunc'] != "unknown" if keep_syn else ~df['ExonicFunc'].isin(
            ["unknown", "synonymous SNV"])
        # (df['AAChange'] != "UNKNOWN") & df['AAChange'].notna()  # keep for splicing
        aa_change = True
        function = ~df['Func'].isin([
            "downstream",
            "intergenic",
            "intronic",
            "ncRNA_exonic",
            "ncRNA_exonic;splicing",
            "ncRNA_intronic",
            "ncRNA_splicing",
            "upstream",
            "upstream;downstream",
            "UTR3",
            "UTR5",
            "UTR5;UTR3"
        ])
        # somatic = df['somatic_status'] != 'Germline'
        return df[exon_func & aa_change & function]

    # filter and sort
    basic_df = filter_basic(anno_df, keep_syn=keep_syn)
    basic_df = sort_df(basic_df)

    # output
    basic_df.to_csv(basic_output, sep='\t', index=False)
    show_output(
        f"Writing basic filtered list ({len(basic_df.index)} muts) to {basic_output}.")

    # ############### FILTER1 ########################################

    def filter1(df, _filter=''):

        # get thresholds from filter_setting_file
        thresh = filter_settings.loc[_filter, :]

        # ##### TUMOR DEPTH ############
        tumor_depth = (df['TR2'] > thresh['variantT']) & (
            df['Tdepth'] > thresh['Tdepth'])

        # ##### VAF ##################
        VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])

        # ##### EB/PoN-Filter ##########
        if thresh['EBscore'] == thresh['EBscore']:
            eb = df['EBscore'] >= thresh['EBscore']
        else:
            eb = True
        pon_eb = (
            eb & (df['PoN-Ratio'] < thresh['PoN-Ratio'])
        ) | (df['PoN-Alt-NonZeros'] < thresh['PoN-Alt-NonZeros'])

        # ##### POPULATION #############
        if thresh['PopFreq'] == thresh['PopFreq']:
            # init a noSNP series with True values for the looping
            noSNP = pd.Series(True, index=df.index)
            # go through pop_cols and filter..
            for col in pop_cols:
                # reformat population columns for filtering
                df.loc[df[col] == ".", col] = 0
                df[col] = df[col].fillna(0).astype(float)
                # combine the looped noSNP with the columns PopFreq checks
                noSNP = noSNP & (df[col] <= thresh['PopFreq'])
        else:
            noSNP = True

        # ## FILTER1 RESCUE ##########
        # per default, rescue all candidate genes
        is_candidate = (df['isCandidate'] == 1) | (
            df['isDriver'] == 1) | (df['ChipFreq'] > 0)

        # ############### AML7 ####################
        # if we are filtering for AML7, we include the 7q genes as interesting
        if filter_name == "filter1-AML7":
            is7q = df['cytoBand'].str.contains('^7q')
            rescue = is7q | is_candidate
        else:
            rescue = is_candidate

        # FINAL FILTER1
        filter_criteria = (tumor_depth & noSNP & pon_eb & VAF) | rescue

        return df[filter_criteria].sort_values(['TVAF'], ascending=False)

    # filter and sort
    filter1_df = filter1(basic_df, _filter='filter1')
    filter1_df = sort_df(filter1_df)

    # write
    show_output(
        f"Writing filter1 list ({len(filter1_df.index)} muts) to {filter1_output}")
    filter1_df.to_csv(filter1_output, sep='\t', index=False)
