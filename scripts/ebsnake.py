import os
import sys

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from ebrun import run_ebscore


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    p = s.params
    c = s.config
    cc = c["EBscore"]

    # load configs into EBconfig
    EBconfig = dict(
        mawk_path=os.path.join(s.scriptdir, "ebscore/shell"),
        genome_split=p.genome_split,
        pon_path=p.pon_path,
        zero_path=os.path.join(p.pon_path, cc["zero_path"]),
        debug=cc["debug"],
        AB_chunk_size=cc["chunksize"]["EBscore"],
        threads=s.threads,
        use_cache=cc["use_cache"],
    )
    # add the EBscore['params'] to EBconfig
    EBconfig.update(cc["params"])

    run_ebscore(
        mut_file=s.input.table,
        tumor_bam=s.input.tumor_bam,
        output_file=str(s.output),
        pon_list=p.pon_list,
        chrom=w.chrom,
        config=EBconfig,
    )


if __name__ == "__main__":
    main(snakemake)