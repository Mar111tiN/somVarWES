import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from script_utils import show_output

# use seaborn plotting defaults
import seaborn as sns

sns.set()


def make_SNP_plot(tumor_file, normal_file):
    # load the tumor file
    t_vaf = pd.read_csv(tumor_file, sep="\t").loc[:, ["FullExonPos", "VAF"]]
    n_vaf = pd.read_csv(normal_file, sep="\t").loc[:, ["FullExonPos", "VAF"]]
    # merge for corresponding SNP pos
    t_n = (
        t_vaf.merge(n_vaf, on="FullExonPos").drop("FullExonPos", axis=1)
        # .query("(VAF_x > 0.05 or VAF_y > 0.05) and (VAF_x < 0.95 or VAF_y < 0.95)")
    )

    sample_t = tumor_file.split("/")[1].replace(".snp", "").replace("_", "")
    sample_n = normal_file.split("/")[1].replace(".snp", "").replace("_", "")
    fig, ax = plt.subplots(figsize=(10, 10))
    _ = ax.scatter(t_n["VAF_x"], t_n["VAF_y"], s=0.2, alpha=0.2)
    _ = ax.set_xlabel(sample_t, fontsize=20)
    _ = ax.set_ylabel(sample_n, fontsize=20)
    # calculate offRate
    df0 = t_n[(t_n > 0.1).any(axis=1)]
    n = len(df0.index)
    df1 = df0[np.abs(df0["VAF_x"] - df0["VAF_y"]) > 0.25]
    m = len(df1.index)
    off_ratio = m / n * 100
    _ = ax.set_title(
        f"{sample_t} vs {sample_n} - offRate {round(off_ratio, 1)}", fontsize=30
    )
    return fig, off_ratio


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    i = s.input
    o = s.output
    p = s.params
    c = s.config
    show_output(f"Plotting SNP plot for {w.sample}")
    fig, _ = make_SNP_plot(i.tumor_snp, i.normal_snp)
    fig.savefig(o.plot)


if __name__ == "__main__":
    main(snakemake)
