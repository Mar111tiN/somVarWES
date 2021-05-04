import os
import pandas as pd
from scipy.stats import fisher_exact as fe
import math
import numpy as np
from multiprocessing import Pool
from script_utils import show_output, show_command

#####################################################################
# ###################### FISHER SCORE ################################


def sort_df(df):
    """
    helper for sorting dfs for chromosomes
    """
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]
    df.loc[:, "Chr"] = pd.Categorical(df["Chr"], chrom_list)
    return df.sort_values(["Chr", "Start"])


def get_fisher_exact(row):
    T1plus = row["TR1+"]
    T1min = row["TR1"] - T1plus
    T2plus = row["TR2+"]
    T2min = row["TR2"] - T2plus
    mat = np.matrix([[T1plus, T2plus], [T1min, T2min]])
    fisher_p = fe(mat)[1]
    if fisher_p:
        return round(-10 * math.log(fisher_p, 10), 1)
    return 5000


def compute_FS(df):
    show_command(f"Computing FisherScore for {len(df.index)} lines")
    df["FisherScore"] = df.apply(get_fisher_exact, axis=1)
    return df


def get_FS_col(df, threads):
    split = np.array_split(df, threads)
    df_chunks = []
    pool = Pool(threads)
    df_chunks = pool.map(compute_FS, split)
    return pd.concat(df_chunks)


def fisher(input, output, threads, log):
    """
    main function for multiprocessor fisher score computation
    takes csv with TR1, TR1+, TR2, TR2+
    """

    show_output(f"Computing FisherScore for file {input[0]}")
    df_table = pd.read_csv(input[0], sep="\t")
    df_fisher = get_FS_col(df_table, threads)

    # reduce to important cols
    # add TR1 and TR2 to make merge unique
    cols = list(df_fisher.columns)[:5] + ["FisherScore", "TR1", "TR2"]
    # select columns and sort
    df_fisher = sort_df(df_fisher[cols])

    # write file to filtered
    df_fisher.to_csv(output[0], sep="\t", index=False)
    show_output(
        f"FisherScore for file {input[0]} written to {output[0]}", color="success"
    )
