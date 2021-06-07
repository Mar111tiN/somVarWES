import os
import sys

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from CNV_raw import TN2CNV
from script_utils_CNV import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    i = s.input
    o = s.output
    p = s.params
    c = s.config
    cc = c["CNV"]

    # load configs into CNVconfig
    CNVconfig = dict(
        mawk_path=os.path.join(s.scriptdir, "myCNV/shell"),
        bed_file=p.bedfile,
        gc_split_path=p.gc_split_path,
        genmap_split_path=p.genmap_split_path,
    )
    # add the CNV params to CNVconfig
    CNVconfig.update(cc)

    show_output(f"Computing coverage and BAF profile for {w.sample} from {i.pileup}", time=True)
    cov_df, _ = TN2CNV(
        clean_TN_pileup_file=i.pileup,
        chrom=w.chrom,
        SNP_output=o.snp,
        config=CNVconfig
    )
    show_output(f"Writing coverage track for {w.sample} to {o.cov}")
    cov_df.to_csv(o.cov, sep="\t", index=False, compression="gzip")
    show_output(f"Finished CNV computation for {w.sample} from {i.pileup}", time=True, color="success")

if __name__ == "__main__":
    main(snakemake)
