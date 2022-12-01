import pandas as pd
import numpy as np
import os
import re
import sys
from yaml import CLoader as Loader, load
from script_utils import show_output


def get_refgene(df, config={}):
    """
    reduce the different gene-ref columns to usable data
    collapse to one Func column for filtering
    """

    # RENAMING ##################################
    # # get the ensGene version from first col containing ".ensgene"
    # ensgene_cols = [col for col in df.columns if ".ensGene" in col]
    # if len(ensgene_cols):
    #     ensgene = ensgene_cols[0].split(".")[1]

    # get rename dict for .refGene annotation
    pat = re.compile(r"(.*)\..*Gene.*")

    refgen_dict = {
        col: col.replace(".", "_") for col in df.columns if re.match(pat, col)
    }

    # rename the columns
    df = df.rename(refgen_dict, axis=1)
    # init the collapsed Func and Gene columns
    df.loc[:, "Func"] = ""
    df.loc[:, "Gene"] = ""

    exonic_types = ["exonic", "splicing", "exonic;splicing"]
    if config["keep_UTR"]:
        exonic_types += ["UTR5", "UTR3", "UTR5;UTR3"]

    for func in [col for col in df.columns if col.startswith("Func_")]:
        isExonic = df[func].isin(exonic_types)
        # rename double genes
        gene_col = func.replace("Func", "Gene")
        df.loc[:, gene_col] = df[gene_col].str.replace(r"([-a-zA-Z0-9]+);\1", r"\1", regex=True)
        # fuse Gene cols into one
        df.loc[isExonic, "Gene"] = df["Gene"] + ";" + df[gene_col]

        # only exonic and exonic;splicing have ExonicFunc values
        # transfer these values to Func and all is save
        # get exonicFunc column from Func col
        exonic_func = "Exonic" + func
        # shorten nonsynonymous SNV to nonsynSNV etc
        df.loc[:, exonic_func] = df[exonic_func].str.replace("synonymous ", "syn")

        # Func == "exonic"  --> Func = ExonicFunc
        df.loc[df[func] == "exonic", func] = df[exonic_func]

        # Func=exonic;splicing --> Func = ExonicFunc;splicing
        df.loc[df[func] == "exonic;splicing", func] = df[exonic_func] + ";splicing"

        # drop the func
        df = df.drop(exonic_func, axis=1)
        # fuse Func cols into one
        df.loc[:, "Func"] = df["Func"] + ";" + df[func]

    # merge Func cols into one
    df.loc[:, "Func"] = (
        df["Func"].str.lstrip(";").str.split(";").apply(lambda x: ";".join(set(x)))
    )
    # merge Gene cols into one
    # gene list (gl) into lambda and reduce via set
    df.loc[:, "Gene"] = (
        df["Gene"].str.lstrip(";").str.split(";").apply(lambda gl: ";".join(set(gl)))
    )

    return df


def get_PON_info(df):
    """
    computes values from the PON column
    """

    STRANDSEP = "="
    ADSEP = "<"

    org_cols = list(df.columns)

    # cheat a bad PON for EB = 0 (|1|1|1|1|1|1|1|1|1|1|1|1|1<1|1|1|1|1|1|1|1|1|1|1()
    bad_pon = "<".join(["|".join(["1"] * 30) for i in range(2)])

    df.loc[df["PON+"] == ".", ["PON+", "PON-"]] = [bad_pon, bad_pon]
    df.loc[:, "PONA+"] = df["PON+"].str.split("<").str[0]
    df.loc[:, "POND+"] = df["PON+"].str.split("<").str[1]
    df.loc[:, "PONA-"] = df["PON-"].str.split("<").str[0]
    df.loc[:, "POND-"] = df["PON-"].str.split("<").str[1]
    df.loc[:, "PONA"] = df["PONA+"] + "|" + df["PONA-"]
    df.loc[:, "POND"] = df["POND+"] + "|" + df["POND-"]

    df.loc[:, "PONSum"] = (
        df["POND"].str.split("|").apply(lambda s: np.array(s).astype(int).sum())
    )
    df.loc[:, "PONAltSum"] = (
        df["PONA"].str.split("|").apply(lambda s: np.array(s).astype(int).sum())
    )
    df.loc[:, "PONRatio"] = df.loc[:, "PONAltSum"] / df.loc[:, "PONSum"]
    df.loc[:, "PONAltNonZeros"] = (
        # count 1 as zero in this scenario
        df["PONA"]
        .str.split("|")
        .apply(lambda s: len(s) - s.count("0") - s.count("1"))
    )

    return df.loc[:, org_cols + ["PONAltSum", "PONRatio", "PONAltNonZeros"]]


def apply_clinscore(df, clinscore_file):
    """
    extract, score and realign clinical columns
    """

    # LOAD YAML

    with open(clinscore_file, "r") as stream:
        clinscore_dict = load(stream, Loader=Loader)

    # ##################################################################
    # ############## ICGC29 ###########
    # get icgc version
    icgc = [c for c in df.columns if c.startswith("icgc")][0].split("_")[0]
    show_output(f"Add ICGC derived columns and compute {icgc}_freq", time=False)
    ICGC = (
        df[f"{icgc}_Affected"]
        .str.extract(r"^([0-9]+)/([0-9]+)$")
        .fillna(0)
        .astype("int")
    )
    df[f"{icgc}_freq"] = (ICGC[0] / ICGC[1]).round(4).fillna(".")

    # ##################################################################

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
    show_output(
        f"Add {cosmic90_version} derived columns and compute {cosmic90_version}_freq",
        time=False,
    )

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
    #     "cosmic90_score": 2,  # derived score
    #     "clinvar_score": 2,  # derived score
    #     "icgc29_freq": 2000,  # derived score --> adjust
    # }
    CLINSCORE = clinscore_dict["CLINSCORE"]
    CLINSCORE[f"{cosmic90_version}_score"] = CLINSCORE["cosmic90_score"]
    CLINSCORE.pop("cosmic90_score", None)
    # GET COMBINED ClinScore
    df["ClinScore"] = 0
    show_output("      Combining clinical scores into ClinScore", time=False)
    for col in CLINSCORE.keys():
        show_output(f"            {col}", time=False)
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
        sheet for sheet in sheets if "Driver" not in sheet and "CHIP" not in sheet
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
    if "CHIP" in sheets:
        hotspots = xl.parse("CHIP")
        df = df.merge(
            hotspots.drop("Protein", axis=1),
            how="left",
            on=["Chr", "Start", "Ref", "Alt", "Gene"],
        )

    return df


def addGenmap(df, genmap_path="", modes=["30_0", "50_0", "75_1", "100_2"], chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]):
    """
    adds the genmap data to the coverage data
    selects only the columns that are given in modes
    choose from [
        '30_0', '30_1', '30_2',
        '50_0', '50_1', '50_2',
        '75_0','75_1', '75_2', '75_3',
        '100_0', '100_1', '100_2', '100_4',
        '150_0','150_1', '150_2', '150_4']

    """

    show_output("Adding genmap data.", time=True)
    
    # sorting for categorical chrom should not be needed as the for loop enforces chrom order
    # df["Chr"] = pd.Categorical(df["Chr"], chrom_list)

    # save cols
    cols = df.columns
    # cycle through chroms
    map_dfs = []
    for chrom in chrom_list:
        # load genmap data for that chromosome
        genmap_file = os.path.join(genmap_path, f"hg38_genmap.{chrom}.txt.gz")
        chrom_df = df.query("Chr == @chrom")
        # check
        if not os.path.isfile(genmap_file):
            show_output(f"Could not find genmap file {genmap_file}", color="warning")
            return

        try:
            cols = ["Chr", "Pos"] + modes
            genmap_df = (
                pd.read_csv(genmap_file, sep="\t", compression="gzip")
                .loc[:, cols]
                .fillna(method="ffill").fillna(method="bfill")
            )
            # rename the mappability cols
            mode_rename = {mode: "map" + mode for mode in modes}
            genmap_df = genmap_df.rename(mode_rename, axis=1)
            show_output(
                f"Loading mappability data for {chrom} from {genmap_file}", time=False
            )
        except Exception as e:
            show_output(
                f"{e}: Could not load mappability data for {chrom} from {genmap_file}",
                color="warning",
            )
            return
        # merge the chrom_df with genmap file
        chrom_df = chrom_df.merge(
            genmap_df, left_on=["Chr", "Start"], right_on=["Chr", "Pos"], how="left"
        ).drop("Pos", axis=1)
        map_dfs.append(chrom_df)
        del chrom_df
    map_df = pd.concat(map_dfs)
    del map_dfs
    show_output("Genmap data added.", time=True, color="success")
    return map_df


def interpolate(df, data_col, ref_col="Start", expand_limit=20):
    """
    interpolates missing values in data_col using linear interpolation based on ref_col
    """
    cols = list(df.columns)
    # set FullExonPos as index for the interpolation method to work on proper intervals
    df = df.reset_index(drop=False).set_index(ref_col, drop=False)
    df.loc[:, data_col] = df[data_col].interpolate(
        method="values", limit=expand_limit, limit_direction="both"
    )
    return df.set_index("index")[cols]


def addGCratio(df, gc_path="", mode="100-10", chrom_list=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]):
    """
    adds the GCdata to the coverage data
    gc_path is path to chrom-split gc-data
    """

    show_output("Adding GCratio data")
    # cycle through chroms
    gc_dfs = []
    for chrom in chrom_list:
        chrom_df = df.query("Chr == @chrom")

        # load gc data for that chromosome
        gc_file = os.path.join(gc_path, f"{chrom}.gc{mode}.gz")
        if not os.path.isfile(gc_file):
            show_output(f"Could not find GC file {gc_file}", color="warning")
            return

        try:
            show_output(f"Loading GC data for {chrom} from {gc_file}", time=False)
            gc_df = pd.read_csv(gc_file, sep="\t", compression="gzip").rename(
                {"GC/AT": "GCratio"}, axis=1
            )
            # add the Chr column for the correct merge
            gc_df.loc[:, "Chr"] = chrom
        except Exception as e:
            show_output(
                f"{e}: Could not load GC data for {chrom} from {gc_file}",
                color="warning",
            )
            return
        # gc data has columns Start -> GC/AT
        chrom_df = chrom_df.merge(gc_df, on=["Chr", "Start"], how="outer")
        
        # print("chrom_df", round(sys.getsizeof(chrom_df)/1024**2), "MB")
        del gc_df
        # # interpolate and remove GCdata rows
        chrom_df = interpolate(chrom_df, "GCratio", ref_col="Start").query("End == End")
        gc_dfs.append(chrom_df)
        # print("gc_dfs", sys.getsizeof(gc_dfs))
        del chrom_df
    # concat chrom_dfs
    df = pd.concat(gc_dfs)
    del gc_dfs
    show_output("GCratio data added!", color="success")
    return df


def rearrange_cols(df, file_name="cols", folder=""):
    """
    use an edit file to filter and rearrange the output columns
    """

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
