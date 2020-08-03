from script_utils import run_cmd


def run_annovar(i, o, annovar, annoparams, annoinfo):
    output_base = o[0].replace('.csv', '')
    # remove header because otherwise annovar will use header as first row
    i_nohead = f"{output_base}.nohead.csv"
    run_cmd(f"mawk 'NR > 1 {{print}}' < {i} > {i_nohead}")
    # &>{log}
    anno_cmd = f"{annovar}/table_annovar.pl {i_nohead} --outfile {output_base} {annoparams}"
    run_cmd(anno_cmd)

    format_cmd = f"cat {output_base}.*.txt | {annoinfo} > {o}"
    run_cmd(format_cmd)

    # cleanup
    run_cmd(f"rm -f {i_nohead}")
    # shell(f"rm -f {output_base}.*_multianno.txt")

    print(f"Annotated with annovar and written to {o}")
