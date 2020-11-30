import os
from varscan_core import convert_varscan2table


def main(s):
    config = s.config
    p = s.params

    build = config["ref"]["build"]
    ref = os.path.join(config["paths"]["mystatic"], config["ref"][build]["genome"])

    convert_varscan2table(
        s.input,
        s.output,
        refgen=ref,
        isVCF=config["varscan"]["vcf"],
        vcf2csv=p.vcf2csv,
        editcsv=p.editcsv,
        coords2annovar=p.coords2annovar,
        varscan2table=p.varscan2table,
    )


if __name__ == "__main__":
    main(snakemake)
