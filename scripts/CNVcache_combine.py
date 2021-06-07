import os
import pandas as pd
import sys

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from combineCNV import make_PON_coverage, make_PON_snp
from script_utils_CNV import show_output


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config
    cc = c["CNV"]

    # load configs into EBconfig
    CNVconfig = dict(
        PON_path="."
    )
    # add the EBscore['params'] to EBconfig
    CNVconfig.update(cc)

    # make the PON_files and save
    show_output("Combining PON data for all chroms and performing normalizations", time=True)
    _, _ = make_PON_coverage(config=CNVconfig, save=True)
    _ = make_PON_snp(config=CNVconfig, save=True)
    show_output("Combining PON data finished", color="success", time=True)


if __name__ == "__main__":
    main(snakemake)
