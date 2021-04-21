import os
import sys

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "HDRdetect/code"))
from HDR_run import HDR_master
from script_utils import show_output


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    i = s.input
    o = s.output
    p = s.params
    c = s.config
    cc = c["HDR"]

    # load all HDR params into HDR_config
    HDR_config = dict(
        mawk_path=os.path.join(s.scriptdir, "HDRdetect/shell"),
        MINQ=c["mpileup"]["Q"],
        MINq=c["mpileup"]["MAPQ"],
        genome_split=p.genome_split,
    )

    # add the HDR['params'] to HDRconfig
    HDR_config.update(cc["params"])

    # get the output file
    out_file = o.HDR_table
    bam_file = i.bam

    # ## run the main HDR function
    HDR_df = HDR_master(
        mut_file=i.filter_file,
        bam_file=bam_file,
        chrom=w.chrom,
        threads=cc["threads"],
        HDR_config=HDR_config,
        pileup_file=i.pileup,
    )

    # write to output
    HDR_df.to_csv(out_file, sep="\t", index=False)
    show_output(
        f"HDRdetect of {bam_file} finished. Writing results to {out_file}",
        time=True,
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)
