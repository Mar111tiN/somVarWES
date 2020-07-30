from script_utils import run_cmd, shell


def run_annovar(i, o, annovar, annoparams, annoinfo):
    output_base = o[0].replace('.csv', '')
    # remove header because otherwise annovar will use header as first row
    headerless_input = f"{output_base}.nohead.csv"
    run_cmd(f"mawk 'NR > 1 {{print}}' < {i} > {headerless_input}")
    # &>{log}
    anno_cmd = f"{annovar}/table_annovar.pl {headerless_input} --outfile {output_base} {annoparams}"
    run_cmd(anno_cmd)

    format_cmd = f"cat {output_base}.*.txt | {annoinfo} > {o}"
    run_cmd(format_cmd)

    # cleanup
    shell(f"rm {headerless_input}")
    # shell(f"rm -f {output_base}.*_multianno.txt")

    print(f"Annotated with annovar and written to {o}")
