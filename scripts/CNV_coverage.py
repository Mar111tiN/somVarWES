import os

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from coverage import get_coverage
from script_utils_CNV import show_output


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    output = o.bedCov

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    CNVconfig = dict(
        mawk_path=os.path.join(s.scriptdir, "myCNV/shell"), bedfile=p.bedfile
    )

    CNVconfig.update(c["CNV"]["coverage"])

    # get the bam_file from the pon_list
    show_output(f"Calculating coverage of {w.sample} on chrom {w.chrom}!", time=True)
    # run the coverage tool
    cov_df = get_coverage(
        i.bam,
        chrom=w.chrom,
        config=CNVconfig,
    )

    show_output(
        f"Writing coverage of {w.sample} on chrom {w.chrom} to {output}",
        color="success",
    )
    cov_df.to_csv(output, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main(snakemake)
