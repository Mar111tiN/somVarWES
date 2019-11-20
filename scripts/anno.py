from os import system as shell
import pandas as pd

annovar = snakemake.config['tools']['annovar']

w = snakemake.wildcards
file_name = f"{w.sample}_{w.tumor_norm}"
input = snakemake.input
output = snakemake.output
params = snakemake.params
anno_params = params.anno
anno_info = params.anno_info
log = snakemake.log

output_base = output[0].replace('.csv', '')
# remove header
headerless_input = f"{output_base}.nohead.csv"
shell(f"mawk 'NR > 1 {{print}}' < {input} > {headerless_input}")
anno_cmd = f"{annovar}/table_annovar.pl {headerless_input} --outfile {output_base} {anno_params}"  # &>{log}
print(anno_cmd)
shell(anno_cmd)
shell(f"rm {headerless_input}")


format_cmd = f"cat {output_base}.*.txt | {anno_info} > {output}"
print(format_cmd)
shell(format_cmd)
shell(f"rm -f {output_base}.*.multianno.txt")

print(f"Annotated with annovar and written to {output}")   
