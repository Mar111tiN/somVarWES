import pandas as pd
import os
from script_utils import show_output, sort_df

sort_cols = {"Chr": True, "TVAF": False, "NVAF": False, "Start": True}


# ######## FILTER 2 ######################
def filter2(df, stringency="", thresh={}, config={}):

    # get thresholds
    show_output(f"filter: filter2-{stringency}")
    # remove syngeneic mutations if keep_syn is active (only applied to filter1
    # .. must be removed now)
    # set to True pd.Series with cheat

    if config["keep_syn"]:
        is_exonic = ~filter1_df["ExonicFunc"].isin(["unknown", "synonymous SNV"])
        filter1_df = filter1_df[is_exonic]
        isSynonymous = df["Start"] > 0
        # go through all the refGene "ExonicFunc" cols (including ensgene)
        for exonic_col in [col for col in df.columns if col.startswith("ExonicFunc")]:
            # set isSynonymous to TRUE if ALL refgene say "synonymous SNV"
            isSynonymous = isSynonymous & (df[exonic_col] == "synonymous SNV")
        df = df.loc[~isSynonymous]

    # DEFINE CANDIDATE
    # used for rescue and special thresholds
    is_candidate = (
        (df["isCandidate"] == 1) | (df["AMLDriver"] == 1) | (df["ChipFreq"] > 0)
    )

    # #### SWITCH FOR AML7
    if "AML7" in config["filter_name"]:
        show_output("Including mutations on 7q to candidate list.", time=False)
        is7q = df["cytoBand"].str.contains("^7q")
        is_candidate = is_candidate | is7q

    # ##### TUMOR DEPTH ############
    tumor_depth = (df["TR2"] >= thresh["variantT"]) & (df["Tdepth"] >= thresh["Tdepth"])

    # ##### VAF ##################
    # #### TVAF
    # either cand with lower TVAF or threshold TVAF
    TVAF = (is_candidate & (df["TVAF"] >= thresh["TVAF4Cand"])) | (
        df["TVAF"] >= thresh["TVAF"]
    )
    # ##### NVAF
    # NVAF is computed from upper threshold and a max proximity to TVAF (VAFSim)
    NVAF = (df["NVAF"] <= (thresh["VAFSim"] * df["TVAF"])) & (
        df["NVAF"] <= thresh["NVAF"]
    )

    # ##### EB/PoN-Filter ##########
    eb = (df["EBscore"] >= thresh["EBscore"]) if thresh["EBscore"] else True

    pon_eb = (eb & (df["PONRatio"] < thresh["PONRatio"])) | (
        df["PONAltNonZeros"] < thresh["PONAltNonZeros"]
    )

    # ############ HDR ####################
    HDR = (df["TumorHDRcount"] <= thresh["HDRcount"]) & (
        df["NormalHDRcount"] <= thresh["HDRcount"]
    )

    # ##### POPULATION #############
    if thresh["PopFreq"] == thresh["PopFreq"]:
        # init a noSNP series with True values for the looping
        noSNP = pd.Series(True, index=df.index)
        # go through pop_cols and filter..
        for col in config["pop_cols"]:
            # reformat population columns for filtering
            df.loc[df[col] == ".", col] = 0
            df[col] = df[col].fillna(0).astype(float)
            # combine the looped noSNP with the columns PopFreq checks
            noSNP = noSNP & (df[col] <= thresh["PopFreq"])
    else:
        noSNP = True

    # ####### STRANDBIAS / POLARITY ##########################
    # Strand Ratio (as FisherScore and simple)
    no_strand_bias = df["FisherScore"] <= thresh["FisherScore"]
    # Strand Polarity (filters out very uneven strand distribution of alt calls)
    if thresh.get("strandPolarity", None):
        pol = thresh["strandPolarity"]
        no_strand_polarity = no_strand_polarity = (df["TR2+"] / df["TR2"] <= pol) & (
            df["TR2+"] / df["TR2"] >= (1 - pol)
        )
    else:
        no_strand_polarity = True

    # let's check whether and is not too hard
    strandOK = no_strand_bias & no_strand_polarity

    # ########## RESCUE #####################
    # Clin_score is used for rescue of all mutations
    clin_score = df["ClinScore"] >= thresh["ClinScore"]
    rescue = clin_score

    # ########### COMBINE CRITERIA ###############
    # rescue is only applied to disputable values within parens
    # criteria outside of parens are hard-filtered
    filter_criteria = (
        tumor_depth & pon_eb & NVAF & (noSNP & strandOK & TVAF & HDR | rescue)
    )

    filter2_df = df[filter_criteria]
    show_output(f"{stringency} {len(filter2_df.index)}", time=False)
    dropped_candidates_df = df[~filter_criteria & is_candidate]
    list_len = len(filter2_df.index)

    return (
        sort_df(filter2_df, cols=sort_cols),
        sort_df(dropped_candidates_df, cols=sort_cols),
        list_len,
    )


def get_filter2(mut_file, filter2_output, filterbam_output=None, config={}):

    # ########### LOADING FILTERS
    show_output(f"Running '{config['filter_name']}'")

    filter_file = config["filter_file"]
    show_output(f"Loading filter file {filter_file}\t", time=False, end="")
    if "xls" in os.path.splitext(filter_file)[1]:
        filter_settings = pd.read_excel(
            filter_file, sheet_name=config["filter_sheet"], index_col=0
        )[:4]
    else:
        filter_settings = pd.read_csv(filter_file, sep="\t", index_col=0)
    show_output("Done.", time=False)

    show_output(f"Loading mutation file {mut_file}. ", time=False, end="")
    filter1_df = pd.read_csv(mut_file, sep="\t", low_memory=False)
    show_output(f"{len(filter1_df.index)} mutations found.", time=False)

    output_base = filter2_output.replace(".loose.csv", "")
    # ################ OUTPUT #############################################################

    filter2_dfs = {}
    dropped_dfs = {}
    df_lengths = {}
    for stringency in ["loose", "moderate", "strict"]:
        (
            filter2_dfs[stringency],
            dropped_dfs[stringency],
            df_lengths[stringency],
        ) = filter2(
            filter1_df,
            stringency=stringency,
            thresh=filter_settings.loc[f"filter2-{stringency}"],
            config=config,
        )
        output_file = f"{output_base}.{stringency}.csv"
        show_output(
            f"Writing filter2.{stringency} ({df_lengths[stringency]}) to {output_file}"
        )
        filter2_dfs[stringency].to_csv(output_file, sep="\t", index=False)
    # write dropped files
    drop_file = f"{output_base}.dropped.csv"
    show_output(
        f"Writing {len(dropped_dfs['loose'].index)} muts to {drop_file}.", time=False
    )
    dropped_dfs["loose"].to_csv(drop_file, sep="\t", index=False)

    if config["excel_output"]:
        excel_file = f"{output_base}.xlsx"
        with pd.ExcelWriter(excel_file) as writer:
            # filter1
            filter1_df.to_excel(writer, sheet_name="filter1", index=False)
            for stringency in ["loose", "moderate", "strict"]:
                filter2_dfs[stringency].to_excel(
                    writer, sheet_name=stringency, index=False
                )

            show_output(f"Writing combined filters to excel file {excel_file}.")
            # write dropped files
            dropped_dfs["loose"].to_excel(writer, sheet_name="dropped", index=False)

    # create the filterbam_table for the selected stringency to be used by filterbam
    if filterbam_output:
        if config["filterbam_stringency"] in filter2_dfs:
            filter2_dfs[config["filterbam_stringency"]].to_csv(
                filterbam_output, sep="\t", index=False
            )
        # if stringency is all or any other, use combination of loose and dropped for filterbam
        else:
            all_df = pd.concat(
                [filter2_dfs["loose"], dropped_dfs["loose"]]
            ).drop_duplicates()
            all_df.to_csv(filterbam_output, sep="\t", index=False)
