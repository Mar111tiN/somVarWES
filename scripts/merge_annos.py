import pandas as pd
from script_utils import sort_df
from derive_cols import rearrange_cols

# ############## SNAKEMAKE ################################


def main(s):
    config = s.config

    i = s.input
    p = s.params
    anno = i.annovar
    fisher = i.fisher
    eb = i.eb_score
    output = str(s.output)
    info_folder = p.info_folder

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

    # fisher_cols = list(fisher_df.columns[5:6])
    # new_cols += fisher_cols  # -->COLS

    if config["EBscore"]["run"]:
        print(f"Loading EBScore from file {eb} and merging into annotation.")
        eb_df = (
            pd.read_csv(eb, sep="\t", dtype={"Chr": str, "Start": int, "End": int})
            .fillna(".")
            .sort_values(["Chr", "Start"])
        )
        # eb_cols = list(eb_df.columns[5:])
        # avoid duplication of duplicate varscan calls
        anno_df = (
            anno_df.merge(eb_df, axis=1), on=["Chr", "Start", "End", "Ref", "Alt"])
            .drop_duplicates()
            .rename({"EB": "EBscore"}, axis=1)
        )

    # ############## REARRANGE COLUMNS  and sort rows ###################################
    anno_df = rearrange_cols(anno_df, file_name="raw_cols", folder=info_folder)
    anno_df = sort_df(anno_df)

    # ############# WRITE OUTPUT ########################################################
    anno_df.to_csv(output, sep="\t", index=False)
    print(f"Writing combined annotation to {output}.")


if __name__ == "__main__":
    main(snakemake)
