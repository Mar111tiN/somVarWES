import os
from varscan_core import convert_varscan2table


def main(s):
    cc = s.config["varscan"]
    p = s.params

    # create config for convert_function
    config = {
        "refgen": p.refgen,
        "isVCF": cc["vcf"],
        "mawk_path": os.path.join(s.scriptdir, "shell"),
    }

    convert_varscan2table(s.input, s.output, config=config)


if __name__ == "__main__":
    main(snakemake)
