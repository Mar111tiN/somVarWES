import os
import sys
# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "ebscore/code"))
from file2matrix import PON2matrix


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
        bed_file=p.bedfile,
        pon_path=p.pon_path,
        temp_dir=os.path.join(p.pon_path, "temp"),
        use_cache=cc["use_cache"],
    )
    # add the EBscore['params'] to EBconfig
    EBconfig.update(cc["params"])

    PON2matrix(
        pon_list=p.pon_list,
        chrom=w.chrom,
        config=EBconfig,
    )


if __name__ == "__main__":
    main(snakemake)