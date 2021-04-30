import os
import sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from ebcache import PONmatrix2AB_multi
from script_utils import show_output


def save_split_df(basename, split_tuple):
    """
    receives (split, split_df) tuple and does the writing
    """
    outfile = basename.replace(".splitcount.", f".{split_tuple[0]}.")
    split_tuple[1].to_csv(outfile, sep="\t", index=False, compression="gzip")
    show_output(f"Written split PON matrix to file {outfile}", multi=True)


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    p = s.params
    c = s.config

    input_file = str(s.input.matrix)
    # load configs into EBconfig
    split = p.split
    base_name = s.output[0].replace(".0.", ".splitcount.")

    # loading the matrix cache file
    show_output(f"Loading cache matrix file for splitting into {p.split} parts.")
    pon_matrix_df = pd.read_csv(
        input_file,
        sep="\t",
        compression="gzip",
    )

    # get the split pon_matrix_df for that split
    split_matrix_df = np.array_split(pon_matrix_df, p.split)
    del pon_matrix_df

    # write data in multicore fashion
    split_pool = Pool(s.threads)
    split_pool.map(
        partial(save_split_df, base_name),
        enumerate(split_matrix_df),
    )
    split_pool.close()
    split_pool.join()

    show_output(
        f"Finished splitting cache matrix file for {w.chrom} into {p.split} parts.",
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)
