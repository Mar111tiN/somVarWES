import os
import sys

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from script_utils_CNV import show_output
from combineCNV import combine_sample_CNV, filter_cov
from rollingCov import rolling_coverage
from rollingSNP import rolling_snp, remergeCNV
from plot import plot_snp, plot_cov, plot_CNV
from plot_params import (
    covfig_params,
    snpfig_params,
    cnvfig_params,
    vaf1,
    vaf2,
    log1,
    log2,
    L2Rmean1,
    L2Rmean2,
)


def main(s):
    """
    snakemake wrapper for CNV combine and rolling metrices
    """

    w = s.wildcards
    o = s.output
    p = s.params
    c = s.config
    cc = c["CNV"]

    sample = f"{w.sample}_{w.tumor}-{w.normal}"

    CNVconfig = dict(
        cov_path="cnv",  # path containing cov.gz files for this sample
        snp_path="cnv",  # path containing snp files for this sample
        PON_path=p.pon_path,
        chrom_list=p.chrom_list
    )

    # add the CNV params to CNVconfigc
    CNVconfig.update(cc)

    # run the function
    show_output(f"Combining CNV data for sample {sample}", time=True)
    cov_df, snp_df = combine_sample_CNV(sample, config=CNVconfig)

    show_output(
        f"Finished combining chrom data for sample {sample}. Writing to files..",
        time=True,
        color="success",
    )

    # ########## raw figures ###################
    cov_df.to_csv(o.cov, sep="\t", index=False)
    snp_df.to_csv(o.snp, sep="\t", index=False)

    # # ########## raw figures ###################
    plot_base = o.plot.replace(f"_{w.tumor}.", "_type.")
    fig, _, _, _ = plot_snp(
        snp_df, plots=[vaf2], chroms="all", region="", **snpfig_params
    )
    fig.savefig(plot_base.replace("type", w.tumor))

    fig, _, _, _ = plot_snp(
        snp_df, plots=[vaf1], chroms="all", region="", **snpfig_params
    )
    fig.savefig(plot_base.replace("type", w.normal))

    fig, _, _, _ = plot_cov(
        cov_df, plots=[log2], chroms="all", region="", **covfig_params
    )
    plot_base = plot_base.replace("snp", "cov")
    fig.savefig(plot_base.replace("type", w.tumor))

    fig, _, _, _ = plot_cov(
        cov_df, plots=[log1], chroms="all", region="", **covfig_params
    )
    fig.savefig(plot_base.replace("type", w.normal))

    show_output(
        f"Finished writing output for sample {sample}.", time=True, color="success"
    )

    # ########## rolling coverage ###############
    show_output(f"Filtering and rolling coverage data for sample {sample}", time=True)
    cov_df = filter_cov(cov_df, config=CNVconfig)
    cov_df = rolling_coverage(cov_df, config=CNVconfig)
    show_output(
        f"Finished rolling coverage for sample {sample}", time=True, color="success"
    )

    show_output(f"Filtering and rolling SNP data for sample {sample}", time=True)
    snp_df = rolling_snp(snp_df, cov_df, config=CNVconfig)
    show_output(f"Writing CNV data for sample {sample}to file {o.CNV}..")
    snp_df.to_csv(o.CNV, index=False, sep="\t", compression="gzip")
    show_output(f"Done!")
    show_output(
        f"Remerging CNV data for sample {sample} with full coverage track and writing to {o.fullCNV}",
        time=True,
    )
    cnv_df = remergeCNV(snp_df, cov_df)
    del cov_df

    cnv_df.to_csv(o.fullCNV, index=False, sep="\t", compression="gzip")
    show_output("Done! Saving plots..")

    # ########## CNV figures ###################
    plot_base = plot_base.replace(".cov.raw", ".cnv")

    fig, _, _, _ = plot_CNV(
        cnv_df,
        cov_plots=[log2, L2Rmean2],
        snp_plots=[vaf2],
        chroms="all",
        region="",
        **cnvfig_params,
    )
    CNVfigT_file = plot_base.replace("type", w.tumor)
    fig.savefig(CNVfigT_file)
    fig, _, _, _ = plot_CNV(
        cnv_df,
        cov_plots=[log1, L2Rmean1],
        snp_plots=[vaf1],
        chroms="all",
        region="",
        **cnvfig_params,
    )
    CNVfigN_file = plot_base.replace("type", w.normal)
    fig.savefig(CNVfigN_file)
    del cnv_df

    plot_base = plot_base.replace(".cnv", ".snp.cnv")
    fig, _, _, _ = plot_CNV(
        snp_df,
        cov_plots=[log2, L2Rmean2],
        snp_plots=[vaf2],
        chroms="all",
        region="",
        **cnvfig_params,
    )
    cnvfigT_file = plot_base.replace("type", w.tumor)
    fig.savefig(cnvfigT_file)
    fig, _, _, _ = plot_CNV(
        snp_df,
        cov_plots=[log1, L2Rmean1],
        snp_plots=[vaf1],
        chroms="all",
        region="",
        **cnvfig_params,
    )
    cnvfigN_file = plot_base.replace("type", w.normal)
    fig.savefig(cnvfigN_file)

    show_output(
        f"Combining of coverage and SNP tracks for {sample} into CNV data finished!",
        time=True,
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)
