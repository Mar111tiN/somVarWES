import os
import pandas as pd
from script_utils import show_output, sort_df


#  ############## BASIC FILTER ####################################
def filter_basic(df, config={}):
    """
    basic cutoff based on gene function in Func col
    """

    exon_func = ["deletion", "insertion", "splicing", "stop", "nonsynSNV"]
    # add allowed Funcs UTR if enabled
    if config["keep_UTR"]:
        exon_func += ["UTR"]
    if config["keep_syn"]:
        exon_func += ["synSNV"]

    df = df.loc[
        df["Func"].str.replace("_splicing", "").str.contains("|".join(exon_func)), :
    ]

    # basic cutoff for WES-data
    df = df.loc[df["Tdepth"] > 10, :]
    return df


def filter1(df, thresh={}, config={}):

    # ##### TUMOR DEPTH ############
    tumor_depth = (df["TR2"] > thresh["variantT"]) & (df["Tdepth"] > thresh["Tdepth"])

    # ##### VAF ##################
    VAF = (df["NVAF"] <= thresh["NVAF"]) & (df["TVAF"] >= thresh["TVAF"])

    # ##### EB/PoN-Filter ##########
    if thresh["EBscore"] == thresh["EBscore"]:
        eb = df["EBscore"] >= thresh["EBscore"]
    else:
        eb = True
    pon_eb = (eb & (df["PONRatio"] < thresh["PONRatio"])) | (
        df["PONAltNonZeros"] < thresh["PONAltNonZeros"]
    )

    # ##### POPULATION #############
    if thresh["PopFreq"] == thresh["PopFreq"]:
        # init a noSNP series with True values for the looping
        noSNP = pd.Series(True, index=df.index)
        # go through pop_cols and filter..
        for col in config["pop_cols"]:
            # reformat population columns for filtering
            df.loc[df[col] == ".", col] = 0
            df.loc[:, col] = df[col].fillna(0).astype(float)
            # combine the looped noSNP with the columns PopFreq checks
            noSNP = noSNP & (df[col] <= thresh["PopFreq"])
    else:
        noSNP = True

    # ## FILTER1 RESCUE ##########
    # per default, rescue all candidate genes
    is_candidate = (
        (df["isCandidate"] == 1) | (df["AMLDriver"] == 1) | (df["ChipFreq"] > 0)
    )

    # ############### AML7 ####################
    # if we are filtering for AML7, we include the 7q genes as interesting
    if "AML7" in config["filter_name"]:
        is7q = df["cytoBand"].str.contains("^7q")
        rescue = is7q | is_candidate
    else:
        rescue = is_candidate

    # FINAL FILTER1
    filter_criteria = (tumor_depth & noSNP & pon_eb & VAF) | rescue

    return df.loc[filter_criteria, :].sort_values(["TVAF"], ascending=False)


def filter_basic_filter1(mut_file, basic_output, filter1_output, config={}):

    show_output(f'Running {config["filter_name"]}')
    show_output(f"Loading mutation file {mut_file}.")
    anno_df = pd.read_csv(mut_file, sep="\t")

    filter_file = config["filter_file"]
    show_output(f"Loading filter file {filter_file}")
    if "xls" in os.path.splitext(filter_file)[1]:
        filter_settings = pd.read_excel(
            filter_file, sheet_name=config["filter_sheet"], index_col=0
        )[:4]
    else:
        filter_settings = pd.read_csv(filter_file, sep="\t", index_col=0)

    # use these population columns for checking PopFreq

    show_output(f"    keep_syn= {config['keep_syn']}")
    show_output(f"    keep_UTR= {config['keep_UTR']}")
    show_output(f"Started editing and basic filtering for {mut_file}.")

    # filter and sort
    basic_df = sort_df(filter_basic(anno_df, config=config))

    # output
    basic_df.to_csv(basic_output, sep="\t", index=False)
    show_output(
        f"Writing basic filtered list ({len(basic_df.index)} muts) to {basic_output}."
    )

    # ############### FILTER1 ########################################
    # get thresholds from filter_setting_file
    thresh = filter_settings.loc["filter1"]
    # filter and sort
    filter1_df = sort_df(filter1(basic_df, thresh=thresh, config=config))

    # write
    show_output(
        f"Writing filter1 list ({len(filter1_df.index)} muts) to {filter1_output}"
    )
    filter1_df.to_csv(filter1_output, sep="\t", index=False)
