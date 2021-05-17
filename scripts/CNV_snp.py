import os

# add myCNV package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from heteroSNP import get_heteroSNP
from script_utils_CNV import show_output

# get the run_shell function to be passed to running code
# the snakemake object has to be passed to retrieve the proper scripts folder
# run_shell = set_path('codeCNV', snakemake)


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config

    p = s.params
    w = s.wildcards
    i = s.input
    o = s.output

    output = o.snp

    ########## CONFIG #######################
    # squeeze out the config for get_coverage
    CNVconfig = dict(
        mawk_path=os.path.join(s.scriptdir, "myCNV/shell"),
        bedfile=p.bedfile,
        genome_split_path=p.genome_split,
        SNPdb_path=p.SNPdb_path,
    )

    CNVconfig.update(c["CNV"]["hetSNP"])
    # input files
    bam_file = i.bam
    show_output(f"Detecting SNPs of {w.sample} on chrom {w.chrom}!", time=True)
    # run the coverage tool
    snp_df = get_heteroSNP(
        bam_file,
        chrom=w.chrom,
        config=CNVconfig,
    )
    show_output(
        f"Writing heteroSNP of {w.sample} on chrom {w.chrom} to {output}",
        color="success",
    )
    snp_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main(snakemake)
