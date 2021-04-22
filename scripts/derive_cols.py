import yaml
import re
import pandas as pd
import numpy as np
import os
import openpyxl


def get_PON_info(df):
    """
    computes values from the PON column
    """
    org_cols = list(df.columns)

    df.loc[:, ["PON+"]] = df["PON"].str.split("-").str[0]
    df.loc[:, ["PON-"]] = df["PON"].str.split("-").str[1]
    df.loc[:, ["PONA+"]] = df["PON+"].str.split("=").str[0]
    df.loc[:, ["POND+"]] = df["PON+"].str.split("=").str[1]
    df.loc[:, ["PONA-"]] = df["PON-"].str.split("=").str[0]
    df.loc[:, ["POND-"]] = df["PON-"].str.split("=").str[1]
    df.loc[:, ["PONA"]] = df["PONA+"] + "|" + df["PONA-"]
    df.loc[:, ["POND"]] = df["POND+"] + "|" + df["POND-"]
    df.loc[:, "PONSum"] = (
        df["POND"].str.split("|").apply(lambda s: np.array(s).astype(int).sum())
    )
    df.loc[:, "PONAltSum"] = (
        df["PONA"].str.split("|").apply(lambda s: np.array(s).astype(int).sum())
    )
    df.loc[:, "PONRatio"] = df.loc[:, "PONAltSum"] / df.loc[:, "PONSum"]
    df.loc[:, "PONAltNonZeros"] = (
        df["PONA"]
        .str.split("|")
        .apply(
            lambda s: np.unique(np.array(s).astype(int), return_counts=True)[1][
                1:
            ].sum()
        )
    )
    return df.loc[:, org_cols + ["PONAltSum", "PONRatio", "PONAltNonZeros"]]


def apply_clinscore(df, clinscore_file):
    """
    extract, score and realign clinical columns
    """

    # LOAD YAML

    with open(clinscore_file, "r") as stream:
        clinscore_dict = yaml.load(stream)

    # ##################################################################
    # ############## ICGC29 ###########
    # get icgc version
    icgc = [c for c in df.columns if c.startswith("icgc")][0].split("_")[0]
    print(f"Add ICGC derived columns and compute {icgc}_freq")
    ICGC = (
        df[f"{icgc}_Affected"]
        .str.extract(r"^([0-9]+)/([0-9]+)$")
        .fillna(0)
        .astype("int")
    )
    df[f"{icgc}_freq"] = (ICGC[0] / ICGC[1]).round(4).fillna(".")

    # ##################################################################
    # ############## COSMIC70 #########
    # get the score dict
    cosmic70_location = clinscore_dict["cosmic70"]["location"]
    print("Add Cosmic70 derived columns and compute cosmic70_score")
    cosmic70_pattern = r"(?:ID=(?P<cosmID>COSM[0-9]+(?:,COSM[0-9]+)?);OCCURENCE=)?(?P<freq>[0-9]+)\((?P<organ>[A-Z_a-z]+)\)"
    df[["cosmic70_ID", "cosmic70_freq", "cosmic70_type"]] = (
        df["cosmic70"]
        .str.extractall(cosmic70_pattern)
        .astype({"cosmID": "str", "freq": "int", "organ": "str"})
        .reset_index("match")
        .drop(columns="match")
        .reset_index()
        .groupby("index")
        .aggregate(
            {
                "cosmID": "min",
                "freq": "sum",
                "organ": lambda col: col.str.cat(sep="+"),
            }
        )
    )
    df.loc[:, "cosmic70_freq"] = df["cosmic70_freq"].fillna(0).astype("int")
    df.loc[:, "cosmic70_ID"] = df["cosmic70_ID"].fillna(".")
    df.loc[:, "cosmic70_type"] = df["cosmic70_type"].fillna(".")
    # drop org column automatically
    df = df.drop(columns="cosmic70")

    def cosmic70score(row):
        """
        row-wise computation of cosmic70 scores
        """
        if row["cosmic70_type"] != ".":
            score = 1
            for location in cosmic70_location.keys():
                if location in row["cosmic70_type"]:
                    score += cosmic70_location[location]
            return row["cosmic70_freq"] * score
        else:
            return 0

    df["cosmic70_score"] = df.apply(cosmic70score, axis=1)

    # ##################################################################
    # ############## COSMIC90-9? #########
    # get the score dict
    cosmic90_type_score = clinscore_dict["cosmic90"]["type"]
    cosmic90_location_score = clinscore_dict["cosmic90"]["location"]
    cosmic90_pattern = (
        r"(?P<count>[0-9]+)x\((?P<types>[^0-9@)]+)@(?P<location>[^0-9@)]+)\)"
    )
    # get cosmic9? version
    cosmic90_version = [c for c in df.columns if c.startswith("cosmic9")][0].split("_")[
        0
    ]
    print(f"Add {cosmic90_version} derived columns and compute {cosmic90_version}_freq")

    def cosmic90_score(row):
        """
        row-wise computation of cosmic90 scores
        """

        return (
            (
                1
                + cosmic90_type_score.get(row["types"], 0)
                + cosmic90_location_score.get(row["location"], 0)
            )
            * (row["types"] != ".")
            * int(row["count"])
        )

    df[f"{cosmic90_version}_score"] = (
        df[f"{cosmic90_version}_type"]
        .str.extractall(cosmic90_pattern)
        .apply(cosmic90_score, axis=1)
        .reset_index()
        .drop(columns="match")
        .groupby("level_0")
        .sum()
        .fillna(0)
    )

    # ##################################################################
    # ############## CLINVAR #################
    CLNDN_factorial = clinscore_dict["clinvar"]["factorial"]
    CLNSIG_score = clinscore_dict["clinvar"]["significance"]

    def get_CLINVARscore(row):
        """
        converts the CLINVAR info into scalar clinvar_score
        """

        def CLNDN2score(clndn):
            """
            accumulates a factor for multiplication with CLNSIG_score
            """

            factor = 1
            for key in CLNDN_factorial:
                if key in clndn:
                    factor *= CLNDN_factorial[key]
            return factor

        if row["CLNDN"] == ".":
            return 0
        return (
            1
            + CLNDN2score(row["CLNDN"])
            + CLNSIG_score.get(row["CLNSIG"].split(",")[0], 0)
        )

    df["clinvar_score"] = df.apply(get_CLINVARscore, axis=1)

    # ############## ====> CLINSCORE #################
    # {
    #     "cosmic70_score": 2,  # derived score
    #     "cosmic91_score": 2,  # derived score
    #     "clinvar_score": 2,  # derived score
    #     "icgc29_freq": 2000,  # derived score --> adjust
    # }
    CLINSCORE = clinscore_dict["CLINSCORE"]
    CLINSCORE[f"{cosmic90_version}_score"] = CLINSCORE["cosmic90_score"]
    CLINSCORE.pop("cosmic90_score", None)
    # GET COMBINED ClinScore
    df["ClinScore"] = 0
    print("      Combining clinical scores into ClinScore")
    for col in CLINSCORE.keys():
        print("            ", col)
        df.loc[df[col] != ".", "ClinScore"] += CLINSCORE[col] * df[df[col] != "."][col]
    return df


def get_gene_lists(df, candidate_excel=""):
    """
    parses an excel file for common genes and mutations and populates list
    """
    xl = pd.ExcelFile(candidate_excel, engine="openpyxl")
    sheets = xl.sheet_names

    # get the candidate scores
    candidates = [
        sheet for sheet in sheets if not "Driver" in sheet and not "CHIP" in sheet
    ]
    df.loc[:, "isCandidate"] = 0
    for cand in candidates:
        # get the dataframe
        cand_df = xl.parse(cand)
        col_name = f"{cand}candidate"
        df = df.merge(cand_df.loc[:, ["Gene", "score"]], on="Gene", how="left").rename(
            {"score": col_name}, axis=1
        )
        df.loc[:, col_name] = df[col_name].fillna(0).astype(int)
        df.loc[:, "isCandidate"] = df.loc[:, "isCandidate"] + df.loc[:, col_name]
    drivers = [sheet for sheet in sheets if "Driver" in sheet]

    # get the drivers
    for driver in drivers:
        # get the dataframe
        driver_df = xl.parse(cand)
        # the sheetname itself will be the col_name
        df = df.merge(
            driver_df.loc[:, ["Gene", "score"]], on="Gene", how="left"
        ).rename({"score": driver}, axis=1)
        df.loc[:, driver] = df[driver].fillna(0).astype(int)

    # get the CHIP mutations
    hotspots = xl.parse("CHIP")
    df = df.merge(
        hotspots.drop("Protein", axis=1),
        how="left",
        on=["Chr", "Start", "Ref", "Alt", "Gene"],
    )

    return df


def rearrange_cols(df, file_name="cols", folder=""):
    # write column list to snakepath/info
    col_file = os.path.join(folder, f"{file_name}.txt")
    # make a dataframe from cols and write to file for inspection
    col_df = pd.DataFrame({"cols": list(df.columns)})
    col_df.to_csv(col_file, index=False, header=False)

    # if an edited file list exists, use that for sorting
    new_col_file = os.path.join(folder, f"{file_name}_edit.txt")
    if os.path.isfile(new_col_file):
        try:
            new_cols = list(pd.read_csv(new_col_file, header=None).iloc[:, 0])
            df = df.loc[:, [col for col in new_cols if col in df.columns]]
        except:
            pass
    # if file does not exist, create it
    else:
        col_df.to_csv(new_col_file, index=False, header=False)
    return df
