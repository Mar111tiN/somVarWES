import os
from anno_core import run_annovar

###### anno script chains several commands for removing and reusing headers


def main(s):
    p = s.params

    def get_shell(tool):
        return os.path.join(s.scriptdir, f"shell/{tool}")

    run_annovar(
        s.input,
        s.output,
        annovar=s.config["tools"]["annovar"],
        annoinfo=get_shell("anno_info.mawk"),
        annoparams=p.anno,
    )


if __name__ == "__main__":
    main(snakemake)
