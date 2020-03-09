import os
import pandas as pd

############# SNAKEMAKE ##################

w = snakemake.wildcards
config = snakemake.config
f_config = config['filter']
filter_name = f_config['filter_name']

mut_file = snakemake.input.filter1


output_base = snakemake.output.filter2.replace('.loose.csv', '')
threads = f_config['threads']
keep_syn = f_config['keep_syn']
filter_file = os.path.join(
    config['paths']['filter_settings'],
    f_config['filter_settings']
)

print(f"Running filter2")
print(f'Loading filter1 file {mut_file}.')
filter1_df = pd.read_csv(mut_file, sep='\t')

print(f"Loading filter file {filter_file}")
filter_settings = pd.read_csv(filter_file, sep='\t', index_col=0)


def filter2(df, _filter='moderate'):
    '''
    creates filtered output using the daniel filtering
    input: pd.dataframe of unfiltered annovar output
    output:
    - filtered/sample_tumor-normal_daniel.csv
    '''
    # get thresholds
    print("filter: ", f"filter2-{_filter}")
    thresh = filter_settings.loc[f"filter2-{_filter}", :]
    tumor_depth = (df['TR2'] > thresh['variantT']) & (
        df['Tdepth'] > thresh['Tdepth'])
    if thresh['EBscore']:
        eb = df['EBscore'] >= thresh['EBscore']
    else:
        eb = True

    pon = (df['PoN-Ratio'] < thresh['PoN-Ratio']
           ) | (df['PoN-Alt-NonZeros'] < thresh['PoN-Alt-NonZeros'])

    # minimum TVAF if not candidate
    is_candidate = (df['isCandidate'] == 1) | (df['isDriver'] == 1)
    no_noise = is_candidate | (df['TVAF'] > 0.05)

    # Strand Ratio (as FisherScore and simple)
    no_strand_bias = df['FisherScore'] <= thresh['FisherScore']

    # Strand Polarity (filters out very uneven strand distribution of alt calls)
    if thresh.get(['strand_polarity'], None):
        pol = thresh['strand_polarity']
        no_strand_polarity = (
            df['TR2+'] <= df['TR2'] - pol) & (df['TR2+'] >= pol)
    else:
        no_strand_polarity = True

    strandedness = no_strand_bias & no_strand_polarity

    # VAF is simple
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >= thresh['TVAF'])

    # Clin_score is used for rescue of all mutations
    clin_score = df['ClinScore'] >= thresh['ClinScore']

    # apply filters to dataframe
    df = df[(tumor_depth & (pon | eb) & strandedness & VAF & no_noise)
            | clin_score].sort_values(['TVAF'], ascending=False)
    list_len = len(df.index)
    return df, list_len


################# OUTPUT #############################################################
print(f"Writing filter2 lists to {output_base}.<stringency>.csv")

excel_file = f"{output_base}.xlsx"
with pd.ExcelWriter(excel_file) as writer:
    filter2_dfs = {}
    df_lengths = {}
    for stringency in ['loose', 'moderate', 'strict']:
        filter2_dfs[stringency], df_lengths[stringency] = filter2(
            filter1_df, _filter=stringency)
        print(f"{stringency}: {df_lengths[stringency]}")
        output_file = f"{output_base}.{stringency}.csv"
        filter2_dfs[stringency].to_csv(output_file, sep='\t', index=False)
        filter2_dfs[stringency].to_excel(
            writer, sheet_name=f'{stringency} <{df_lengths[stringency]}>', index=False)
    print(f"Writing combined filters to excel file {excel_file}.")
    filter1_df.to_excel(
        writer, sheet_name=f'filter1 <{len(filter1_df.index)}>', index=False)


# Writing mutation list for filterbam
stringency = config['filter_bam']['stringency_for_bam']
list_for_bam = snakemake.output.filter2_for_filterbam
# select stringency based on filter_bam config
filter2_dfs[stringency].to_csv(list_for_bam, sep='\t', index=False)
