import os
import sys
import pandas as pd
import numpy as np

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from ebcache import PONmatrix2AB_multi
from script_utils import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    p = s.params
    c = s.config
    cc = c["EBscore"]

    output = str(s.output)
    input = str(s.input.matrix)
    # load configs into EBconfig
    EBconfig = dict(
        pon_path=p.pon_path,
        threads=s.threads,
        temp_dir=os.path.join(p.pon_path, "temp"),
        chunk_size=cc["AB_chunk_size"]["EBcache"],
    )
    # add the EBscore['params'] to EBconfig
    EBconfig.update(cc["params"])

    # loading the matrix cache file
    show_output(f"Loading cache matrix file {input}.")
    pon_matrix_df = pd.read_csv(input, sep="\t", compression="gzip")

    split = int(w.split)
    # get the split pon_matrix_df for that split
    split_pon_matrix_df = np.array_split(pon_matrix_df, cc["ABcache_split"])[split]
    del pon_matrix_df
    matrix_len = len(split_pon_matrix_df.index)
    show_output(
        f"Finished. Computing ABcache for {matrix_len} positions ({matrix_len * 12} lines).",
        time=False,
    )
    AB_df = PONmatrix2AB_multi(
        split_pon_matrix_df,
        config=EBconfig,
    )

    AB_df.to_csv(output, sep="\t", compression="gzip")
    show_output(
        f"Finished! Written ABcache split {split} for chrom {w.chrom} to {os.path.join(p.pon_path, output)}",
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)