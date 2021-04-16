import os
from varscan_core import convert_varscan2table
from script_utils import make_mawk


def main(s):
    cc = s.config["varscan"]
    p = s.params

    def get_shell(tool):
        return os.path.join(s.scriptdir, f"shell/{tool}")

    convert_varscan2table(
        s.input, s.output, refgen=p.refgen, isVCF=cc["vcf"], mawk=make_mawk(s)
    )


if __name__ == "__main__":
    main(snakemake)
