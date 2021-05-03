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
    split = w.split
    output = str(s.output)
    input_file = str(s.input.matrix)
    # load configs into EBconfig
    EBconfig = dict(
        zero_path=os.path.join(p.pon_path, cc["zero_path"]),
        threads=s.threads,
        chunksize=cc["chunksize"]["EBcache"],
    )
    # add the EBscore['params'] to EBconfig
    EBconfig.update(cc["params"])

    # loading the matrix cache file
    show_output(f"Loading cache matrix split {input_file} as generator.")
    pon_matrix_gen = pd.read_csv(
        input_file, sep="\t", compression="gzip", chunksize=EBconfig["chunksize"]
    )

    AB_df = PONmatrix2AB_multi(
        pon_matrix_gen,
        config=EBconfig,
    )

    AB_df.to_csv(output, sep="\t", compression="gzip", index=False)
    show_output(
        f"Finished! Written ABcache split {split} for chrom {w.chrom} to {os.path.join(p.pon_path, output)}",
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)
