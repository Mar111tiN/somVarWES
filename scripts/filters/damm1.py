import os
import pandas as pd

############# SNAKEMAKE ##################

w = snakemake.wildcards
config = snakemake.config
f_config = config['filter']
filter_name = f_config['filter_name']
i = snakemake.input[0]
filter1_file = snakemake.output.filter1
filter_basic_file = snakemake.output.basic
threads = f_config['threads']
keep_syn = f_config['keep_syn']
filter_setting_file = os.path.join(
    config['paths']['filter_settings'],
    f_config['filter_settings']
)

# get and run the respective filter as a shell script:
print(f"Running filter {filter_name}")
print(f"    keep_syn= {keep_syn}")

print(f'Started editing and basic filtering for {i}.')
anno_df = pd.read_csv(i, sep='\t')


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


basic_df.to_csv(filter_basic_file, sep='\t', index=False)
print(f"Writing basic filtered list to {filter_basic_file}.")

# ############### FILTER1 ########################################
print(f"Loading filter file {filter_setting_file}")
filter_settings = pd.read_csv(filter_setting_file, sep='\t', index_col=0)

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
    tumor_depth = (df['TR2'] > thresh['variantT']) & (
        df['Tdepth'] > thresh['Tdepth'])

    # EBFilter
    if thresh['EBscore']:
        eb = df['EBscore'] >= thresh['EBscore']
    else:
        eb = True

    # other PanelOfNormal metrices
    pon = (df['PoN-Ratio'] < thresh['PoN-Ratio']
           ) | (df['PoN-Alt-NonZeros'] <= thresh['PoN-Alt-NonZeros'])

    # Strand Bias (as FS)
    no_strand_bias = df['FisherScore'] < thresh['FisherScore']
    VAF = (df['NVAF'] <= thresh['NVAF']) & (df['TVAF'] >=
                                            thresh['TVAF']) & (df['TVAF'] > df['NVAF'])

    return df[tumor_depth & (eb | pon) & no_strand_bias & VAF]


# from Kenichi Data
# misRate_tumor > 0.05 or CoSMIC OCCURRENCE
# depth_tumor > 30 ?
# variantNum_tumor > 3
# EBscore >=4 or CoSMIC OCCURRENCE

filter1_df = filter1(basic_df, _filter='filter1')
print(f"Writing filter1 file to {filter1_file}")
filter1_df.to_csv(filter1_file, sep='\t', index=False)
