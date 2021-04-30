import os
import sys
import pandas as pd

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from zerocache import collapse_zeros
from script_utils import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    p = s.params
    c = s.config
    cc = c["EBscore"]

    # get the pon_size from the PONlist

    pon_file = os.path.join(p.pon_path, p.pon_list)
    ponsize = len(pd.read_csv(pon_file, sep="\t", header=None).index)
    zero_path = os.path.join(p.pon_path, cc["zero_path"])

    for ps in [ponsize, ponsize - 1]:
        collapse_zeros(
            zero_path,
            ponsize=ps,
            ZDfactor=cc["params"]["ZDfactor"],
            reflat=p.reflat,
        )


if __name__ == "__main__":
    main(snakemake)