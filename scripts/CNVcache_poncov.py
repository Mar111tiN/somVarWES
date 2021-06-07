import os
import sys

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from CNV_raw import PON2CNV
from script_utils_CNV import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    o = s.output
    p = s.params
    c = s.config
    cc = c["CNV"]

    # load configs into EBconfig
    CNVconfig = dict(
        mawk_path=os.path.join(s.scriptdir, "myCNV/shell"),
        bed_file=p.bedfile,
        PON_path=".",
        gc_split_path=p.gc_split_path,
        genmap_split_path=p.genmap_split_path,
    )
    # add the EBscore['params'] to EBconfig
    CNVconfig.update(cc)

    cov_df, _ = PON2CNV(
        chrom=w.chrom,
        config=CNVconfig,
    )

    cov_df.to_csv(o.poncov, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main(snakemake)
