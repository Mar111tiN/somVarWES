import os
import sys
import pandas as pd

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from zerocache import collapse_zeros


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    p = s.params
    c = s.config
    cc = c["EBscore"]

    # get the pon_size from the PONlist

    pon_file = os.path.join(p.pon_path, p.pon_list)
    pon_size = len(pd.read_csv(pon_file, sep="\t", header=None).index)

    collapse_zeros(
        os.path.join(p.pon_path, "zero"),
        pon_size=pon_size,
        zero_condense_factor=cc['params']['zero_condense_factor'],
        reflat=p.reflat
    )


if __name__ == "__main__":
    main(snakemake)