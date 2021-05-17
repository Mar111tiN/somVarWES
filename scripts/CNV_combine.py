import os

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from script_utils_CNV import show_output
from combineCNV import combine_CNV


def main(s):
    """
    snakemake wrapper for CNV combine and rolling metrices
    """

    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    snp_list = i.snp
    cov_list = i.cov

    pon_path = p.pon_path
    sample = f"{w.sample}_{w.type}"
    # run the function
    show_output(f"Combining CNV data for sample {sample}", time=True)
    cov_df, snp_df = combine_CNV(
        snp_list,
        cov_list,
        pon_cov_path=os.path.join(pon_path, "chromCov"),
    )

    show_output(
        f"Finished combining chrom data for sample {sample}.",
        time=True,
        color="success",
    )
    cov_df.to_csv(o.cov, sep="\t", index=False)
    snp_df.to_csv(o.snp, sep="\t", index=False)
    show_output(
        f"Finished writing output for sample {sample}.", time=True, color="success"
    )


if __name__ == "__main__":
    main(snakemake)
