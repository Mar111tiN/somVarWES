import pandas as pd
import re
from script_utils import sort_df

# ############## SNAKEMAKE ################################


def main(s):
    config = s.config

    i = s.input
    anno = i.annovar
    fisher = i.fisher
    eb = i.eb_score
    output = s.output
    threads = s.threads
    log = s.log

    # ############## ANNOVAR FILE ########################################
    # ####################################################################

    print(f"Loading annovar file and adjusting columns {anno}.")
    anno_df = (
        pd.read_csv(
            anno,
            sep="\t",
            low_memory=False,
            dtype={"Chr": str, "Start": int, "End": int},
            na_values=["---", "NaN", "NAN"],
        )
        .fillna(".")
        .sort_values(["Chr", "Start"])
    )
    anno_df["Tdepth"] = anno_df["TR1"] + anno_df["TR2"]
    anno_df["TVAF"] = (anno_df["TR2"] / (anno_df["TR1"] + anno_df["TR2"])).round(3)
    # anno_df = anno_df.query('TVAF = TR1 / (TR1+TR2)')
    anno_df["Ndepth"] = anno_df["NR1"] + anno_df["NR2"]
    anno_df["NVAF"] = (anno_df["NR2"] / (anno_df["NR1"] + anno_df["NR2"])).round(3)

    # RENAMING ##################################

    # get rename dict for .refGene annotation
    refgen_dict = {
        col: col.replace(".refGene", "") for col in anno_df.columns if ".refGene" in col
    }
    # get rename dict for .ensGene annotation
    refgen_dict = {
        col: re.sub(".ensGene[0-9]*", "", col)
        for col in anno_df.columns
        if ".ensGene" in col
    }
    # rename the columns
    anno_df = anno_df.rename(columns=refgen_dict)

    # resort the columns
    cols = list(anno_df.columns)
    start_cols = cols[:11]
    quant_cols = cols[-15:]
    anno_cols = cols[11:-15]

    # MERGE Func into ExonicFunc ##################################
    # merge ExonFunc from Func, if ExonFunc not available
    anno_df.loc[anno_df["ExonicFunc"] == ".", "ExonicFunc"] = anno_df["Func"]

    # ############ -->COLS #####################################
    new_cols = start_cols + quant_cols

    # ############## MERGE WITH EB AND FISHER ############################
    #####################################################################

    print(f"Loading FisherScore from file {fisher} and merging into annotation.")
    fisher_df = (
        pd.read_csv(fisher, sep="\t", dtype={"Chr": str, "Start": int, "End": int})
        .fillna(".")
        .sort_values(["Chr", "Start"])
    )
    # merge with TR1 and TR2 to avoid duplication of duplicate varscan calls
    anno_df = anno_df.merge(
        fisher_df, on=["Chr", "Start", "End", "Ref", "Alt", "TR1", "TR2"]
    )
    fisher_cols = list(fisher_df.columns[5:6])
    new_cols += fisher_cols  # -->COLS
    if config["EBFilter"]["run"]:
        print(f"Loading EBScore from file {eb} and merging into annotation.")
        eb_df = (
            pd.read_csv(eb, sep="\t", dtype={"Chr": str, "Start": int, "End": int})
            .fillna(".")
            .sort_values(["Chr", "Start"])
        )
        eb_cols = list(eb_df.columns[5:])
        # avoid duplication of duplicate varscan calls
        anno_df = anno_df.merge(
            eb_df, on=["Chr", "Start", "End", "Ref", "Alt"]
        ).drop_duplicates()
        new_cols += eb_cols  # -->COLS

    # ############## REARRANGE COLUMNS ###################################
    #####################################################################

    new_cols += anno_cols
    anno_df = sort_df(anno_df[new_cols])
    anno_df.to_csv(output[0], sep="\t", index=False)
    print(f"Writing combined annotation to {output[0]}.")


if __name__ == "__main__":
    main(snakemake)
