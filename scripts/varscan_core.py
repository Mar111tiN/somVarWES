import os
from os import system as run
from script_utils import show_output, run_cmd


def convert_varscan2table(input_files, o, config={}):
    """
    converts the vcf output from varscan into readable tab-sep data
    """
    # unwrapper for the mawk tools
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    output_files = []

    if config["isVCF"]:
        for input in input_files:
            # standard output
            # print('input:', input)
            out = input.replace("vcf", "table")
            out_ln = input.replace("vcf", "ln.vcf")
            # sometimes, varscan outputs degenerate ref bases like W leading to errors in bcftools
            # .. here is a dirty fix:
            mawk_cmd = 'BEGIN {OFS="\t"} $4 == "W" {print $1,$2,$3,"N",$5,$6,$7,$8,$9,$10,$11; next;}{print}'
            run_cmd(
                f"cat {input} | mawk '{mawk_cmd}' > {input}.new; mv {input}.new {input}"
            )
            # for indel realignment, files have to be bgzipped and indexed with tabix
            run_cmd(f"bgzip < {input} > {input}.gz")
            run_cmd(f"tabix {input}.gz")
            run_cmd(f"bcftools norm -f {config['refgen']} -o {out_ln} {input}.gz")
            run_cmd(
                f"cat {out_ln} | sed 's/SOMATIC;//' | {mawk('vcf2csv')} | {mawk('editcsvVarscan')} | {mawk('coords2annovar')} > {out}"
            )
            output_files.append(out)
            # cleanup
            run(f"rm {input}.gz; rm {input}.gz.tbi; rm {input}; mv {out_ln} {input}")

    else:
        for input in input_files:
            # print('input:', input)
            out = f"{input}.table"
            output_files.append(out)
            # varscan output has to be converted to avinput file format
            run_cmd(f"{mawk('varscan2table')} < {input} > {out}")

    # CONCAT THE FILES
    # rm first two lines (the headers of both files)
    run_cmd(
        f"cat {' '.join(output_files)} | sort -V -k1,2 | mawk 'NR > 2 {{ print }}' >  {o}"
    )

    run(f"rm -f {' '.join(output_files)}")
    show_output(
        f"Concated {input_files[0]} and {input_files[1]} into {o} for annotation.",
        color="success",
    )
    # shell(f"rm {' '.join(output_files)}")
