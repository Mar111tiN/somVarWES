import os
from varscan_core import convert_varscan2table


def main(s):
    cc = s.config["varscan"]
    p = s.params

    def get_shell(tool):
        return os.path.join(s.scriptdir, f"shell/{tool}")

    convert_varscan2table(
        s.input,
        s.output,
        refgen=p.refgen,
        isVCF=cc["vcf"],
        vcf2csv=get_shell("vcf2csv.mawk"),
        editcsv=get_shell("editcsvVarscan.mawk"),
        coords2annovar=get_shell("coords2annovar.mawk"),
        varscan2table=get_shell("varscan2table.mawk"),
    )


if __name__ == "__main__":
    main(snakemake)
