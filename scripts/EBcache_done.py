import os
import sys
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

    zdf = collapse_zeros(p.zero_path, zero_condense_factor=cc['params']['zero_condense_factor'])


if __name__ == "__main__":
    main(snakemake)