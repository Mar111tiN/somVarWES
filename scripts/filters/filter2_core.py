import pandas as pd
import os
from script_utils import show_output


def get_filter2(mut_file, filter2_output,
                filter_file, filter_sheet, filter_name, keep_syn=False,
                filterbam_output=None, filterbam_stringency='moderate'):

    show_output(f'Loading mutation file {mut_file}.')
    filter1_df = pd.read_csv(mut_file, sep='\t')
    show_output('Done.', time=False)

    # remove syngeneic mutations if keep_syn is active (only valid for filter1)
    if keep_syn:
        is_exonic = ~filter1_df['ExonicFunc'].isin(["unknown", "synonymous SNV"])
        filter1_df = filter1_df[is_exonic]

    # ########### LOADING FILTERS
    show_output(f"Running filter2 {filter_name}")
    show_output(f"\tLoading filter file {filter_file}", time=False)
    if "xls" in os.path.splitext(filter_file)[1]:
        filter_settings = pd.read_excel(
            filter_file, sheet_name=filter_sheet, index_col=0)[:4]
    else:
        filter_settings = pd.read_csv(filter_file, sep='\t', index_col=0)
    show_output('Done.', time=False)

    # use these population columns for checking PopFreq
    # could be refactored into params
    pop_cols = ['gnomAD_exome_ALL', 'esp6500siv2_all', 'dbSNP153_AltFreq']

    output_base = filter2_output.replace('.loose.csv', '')

    # ######## FILTER 2 ######################
    def filter2(df, _filter='moderate'):

        # get thresholds
        show_output("filter: ", f"filter2-{_filter}")
        thresh = filter_settings.loc[f"filter2-{_filter}", :]

        # DEFINE CANDIDATE
        # used for rescue and special thresholds
        is_candidate = (df['isCandidate'] == 1) | (
            df['isDriver'] == 1) | (df['ChipFreq'] > 0)

        # #### SWITCH FOR AML7
        if "AML7" in filter_name:
            is7q = df['cytoBand'].str.contains('^7q')
            is_candidate = is_candidate | is7q

        # ##### TUMOR DEPTH ############
        tumor_depth = (df['TR2'] >= thresh['variantT']) & (
            df['Tdepth'] >= thresh['Tdepth'])

        # ##### VAF ##################
        # #### TVAF
        # either cand with lower TVAF or threshold TVAF
        TVAF = (is_candidate & (df['TVAF'] >= thresh['TVAF4Cand'])) | (
            df['TVAF'] >= thresh['TVAF'])
        # ##### NVAF
        # NVAF is computed from upper threshold and a max proximity to TVAF (VAFSim)
        NVAF = (df['TVAF'] > ((1 + thresh['VAFSim']) * df['NVAF'])
                ) & (df['NVAF'] <= thresh['NVAF'])

        # ##### EB/PoN-Filter ##########
        eb = (df['EBscore'] >= thresh['EBscore']) if thresh['EBscore'] else True

        pon_eb = (eb & (df['PoN-Ratio'] < thresh['PoN-Ratio'])
                ) | (df['PoN-Alt-NonZeros'] < thresh['PoN-Alt-NonZeros'])

        # ############ HDR ####################
        HDR = (df['TumorHDRcount'] <= thresh['HDRcount']) & (
            df['NormalHDRcount'] <= thresh['HDRcount'])

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

        # ####### STRANDBIAS / POLARITY ##########################
        # Strand Ratio (as FisherScore and simple)
        no_strand_bias = df['FisherScore'] <= thresh['FisherScore']
        # Strand Polarity (filters out very uneven strand distribution of alt calls)
        if thresh.get('strandPolarity', None):
            pol = thresh['strandPolarity']
            no_strand_polarity = no_strand_polarity = (
                df['TR2+'] / df['TR2'] <= pol) & (df['TR2+'] / df['TR2'] >= (1-pol))
        else:
            no_strand_polarity = True

        strandOK = no_strand_bias | no_strand_polarity

        # ########## RESCUE #####################
        # Clin_score is used for rescue of all mutations
        clin_score = df['ClinScore'] >= thresh['ClinScore']
        rescue = clin_score

        # ########### COMBINE CRITERIA ###############
        # rescue is only applied to disputable values within parens
        # criteria outside of parens are hard-filtered
        filter_criteria = tumor_depth & pon_eb & NVAF & (
            noSNP & strandOK & TVAF & HDR | rescue)

        filter2_df = df[filter_criteria]
        show_output(stringency, len(filter2_df.index))
        dropped_candidates_df = df[~filter_criteria & is_candidate]
        list_len = len(filter2_df.index)
        return filter2_df, dropped_candidates_df, list_len

    # ################ OUTPUT #############################################################
    excel_file = f"{output_base}.xlsx"

    with pd.ExcelWriter(excel_file) as writer:
        # filter1
        filter1_df.to_excel(
            writer, sheet_name='filter1', index=False)
        filter2_dfs = {}
        dropped_dfs = {}
        df_lengths = {}
        for stringency in ['loose', 'moderate', 'strict']:
            filter2_dfs[stringency], dropped_dfs[stringency], df_lengths[stringency] = filter2(
                filter1_df, _filter=stringency)
            output_file = f"{output_base}.{stringency}.csv"
            show_output(f"Writing filter2.{stringency} ({df_lengths[stringency]}) to {output_file}")
            filter2_dfs[stringency].to_csv(output_file, sep='\t', index=False)
            filter2_dfs[stringency].to_excel(
                writer, sheet_name=stringency, index=False)

        # write dropped files
        drop_file = f"{output_base}.dropped.csv"
        show_output(f"Writing {len(dropped_dfs['loose'].index)} muts to {drop_file}.", time=False)
        dropped_dfs['loose'].to_csv(drop_file, sep='\t', index=False)

        show_output(f"Writing combined filters to excel file {excel_file}.")

        dropped_dfs['loose'].to_excel(writer, sheet_name='dropped', index=False)

    if filterbam_output:
        filter2_dfs[filterbam_stringency].to_csv(filterbam_output, sep='\t', index=False)
