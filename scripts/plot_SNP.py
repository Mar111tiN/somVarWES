import pandas as pd
import os
import sys

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from plot import make_SNP_plot
from script_utils import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    i = s.input
    o = s.output

    sample = f"{w.sample}_{w.tumor}-{w.normal}"
    show_output(f"Plotting SNP plot for {sample} and saving to")
    cnv_df = pd.read_csv(i.CNV, sep="\t", compression="gzip")
    fig, _ = make_SNP_plot(sample, cnv_df)
    fig.savefig(o.plot)


if __name__ == "__main__":
    main(snakemake)
