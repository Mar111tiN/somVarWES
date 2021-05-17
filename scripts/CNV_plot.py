import os
import pandas as pd

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from script_utils_CNV import show_output
from plot import plot_snp2


fig_params = dict(
    figsize=(24, 8),
    colormap="coolwarm_r",
    color_chroms=True,
    ylim=(0, 1),
    cov_offset=0.1,  # how much log2ratio=0 is shifted above SNP-data
    cov_height=0.5,
    label_size=13,
)

log2 = dict(
    title="log2ratio",
    plot_type="scatter",  # ['line', 'scatter']
    data="log2ratio",
    plot_args=dict(linewidth=0.3, color="black", s=0.8, alpha=0.7),
)

log2mean = dict(
    title="log2ratiomean",
    plot_type="line",  # ['line', 'scatter']
    data="log2ratiomean",
    plot_args=dict(linewidth=1, color="yellow", alpha=0.7),
)

vaf = dict(
    title="VAF",
    plot_type="scatter",  # ['line', 'scatter']
    data="VAF",
    plot_args=dict(s=2, color="black", cmap="viridis", alpha=1),
)


def main(s):
    """
    snakemake wrapper for CNV combine and rolling metrices
    """

    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    sample = f"{w.sample}_{w.type}"
    # run the function
    show_output(f"Plotting CNV data for sample {sample}", time=True)

    snp_df = pd.read_csv(i.snp, sep="\t").query("0.05 < VAF < 0.98")

    fig, _, _, _ = plot_snp2(
        snp_df, snp_plots=[vaf], cov_plots=[log2, log2mean], chroms="all", **fig_params
    )

    fig.savefig(o.plot)
    show_output(f"Finished Plotting for sample {sample}.", time=True, color="success")


if __name__ == "__main__":
    main(snakemake)
