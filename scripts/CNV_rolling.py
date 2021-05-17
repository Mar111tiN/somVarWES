import os
import pandas as pd

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from rollingCov import apply_rolling_coverage
from script_utils_CNV import show_output


def main(s):
    """
    snakemake wrapper for CNV combine and rolling metrices
    """

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output
    sample = w.sample

    # run the function
    show_output(
        f"Loading combined coverage and SNP data for sample {sample}", time=True
    )
    cov_df = pd.read_csv(i.cov, sep="\t")
    snp_df = pd.read_csv(i.snp, sep="\t")
    show_output(f"Performing rolling coverage computation for sample {sample}.")
    snpcov_df, rolling_cov_df = apply_rolling_coverage(
        snp_df, cov_df, config=c["CNV"]["combine"]
    )

    show_output(
        f"Rolling coverage computation for sample {sample} finished.",
        time=True,
        color="success",
    )
    rolling_cov_df.to_csv(o.rcov, sep="\t", index=False)
    snpcov_df.to_csv(o.snp, sep="\t", index=False)
    show_output(
        f"Finished writing output for sample {sample}.", time=True, color="success"
    )
    # apply rolling coverage


if __name__ == "__main__":
    main(snakemake)
