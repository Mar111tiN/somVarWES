import os
import pandas as pd

# ############ SNAKEMAKE ##################

w = snakemake.wildcards
config = snakemake.config
f_config = config['filter']

mut_file = snakemake.input[0]

filter_name = f_config['filter1']
filter_file = os.path.join(
    config['paths']['filter_settings'],
    f_config['filter_settings']
)
sheet = f_config['excel_sheet']
threads = f_config['threads']
keep_syn = f_config['keep_syn']

basic_file = snakemake.output.basic
filter1_file = snakemake.output.filter1

print(f"Running filter1 \"{filter_name}\"")
print(f"Loading mutation file {mut_file}.")
anno_df = pd.read_csv(mut_file, sep='\t')

print(f"Loading filter file {filter_file}")
if "xls" in os.path.splitext(filter_file)[1]:
    filter_settings = pd.read_excel(
        filter_file, sheet_name=sheet, index_col=0)[:4]
else:
    filter_settings = pd.read_csv(filter_file, sep='\t', index_col=0)

# use these population columns for checking PopFreq
# could be refactored into params
pop_cols = ['gnomAD_exome_ALL', 'esp6500siv2_all', 'dbSNP153_AltFreq']

print(f"    keep_syn= {keep_syn}")
print(f'Started editing and basic filtering for {mut_file}.')


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


basic_df = filter_basic(anno_df, keep_syn=keep_syn)

basic_df.to_csv(basic_file, sep='\t', index=False)
print(
    f"Writing basic filtered list ({len(basic_df.index)} muts) to {basic_file}.")

# ############### FILTER1 ########################################

# filter1_setting = {
#     'variantT': 2,
#     'Tdepth': 20,
#     'EBscore': 1,
#     'PoN-Ratio': 0.001,
#     'PoN-Alt-NonZeros': 4,
#     'FisherScore': 50,
#     'TVAF': 0.01,
#     'NVAF': 0.3,
# }


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
    pon_eb = (eb & (df['PoN-Ratio'] < thresh['PoN-Ratio'])
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
    is7q = df['cytoBand'].str.contains('^7q')
    is_candidate = (df['isCandidate'] == 1) | (df['isDriver'] == 1) | (df['ChipFreq'] > 0)
    rescue = is7q | is_candidate

    # FINAL FILTER1
    filter_criteria = (tumor_depth & noSNP & pon_eb & VAF) | rescue

    return df[filter_criteria].sort_values(['TVAF'], ascending=False)


filter1_df = filter1(basic_df, _filter='filter1')
print(f"Writing filter1 file to {filter1_file}")
filter1_df.to_csv(filter1_file, sep='\t', index=False)
