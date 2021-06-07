import os
import sys
import pandas as pd

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))

from script_utils_CNV import show_output
from convert import write_ASCAT


def main(s):
    """
    snakemake wrapper for writing ASCAT output
    """
    w = s.wildcards
    i = s.input

    sample = f"{w.sample}_{w.tumor}-{w.normal}"
    show_output(f"Loading {i.CNV} of sample {sample} for ASCAT conversion")
    cnv_df = pd.read_csv(i.CNV, sep="\t", compression="gzip")
    write_ASCAT(cnv_df, sample=sample, outpath=f"CNV/{w.sample}/pre")


if __name__ == "__main__":
    main(snakemake)
