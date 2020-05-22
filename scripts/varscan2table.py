from os import system as shell
import os
from script_utils import show_output, show_command

w = snakemake.wildcards
config = snakemake.config


input_files = snakemake.input
output = snakemake.output
params = snakemake.params
VCF = config['varscan']['vcf']
vcf2csv = params.vcf2csv
editcsv = params.editcsv
coords2annovar = params.coords2annovar
# pseudopipe the tree awk scripts
vcf2table = f"{vcf2csv} | {editcsv} | {coords2annovar}"
varscan2table = params.varscan2table

log = snakemake.log
build = config['ref']['build']
ref = os.path.join(config['paths']['mystatic'], config['ref'][build]['genome'])

def run_cmd(cmd):
    show_command(cmd, multi=False)
    shell(cmd)

output_files = []
if VCF:  
    for input in input_files:
        # standard output
        # print('input:', input)
        out = input.replace('vcf', 'table')
        out_ln = input.replace('vcf', 'ln.vcf')
        # for indel realignment, files have to be bgzipped and indexed with tabix   
        run_cmd(f"bgzip < {input} > {input}.gz")
        run_cmd(f"tabix {input}.gz")
        run_cmd(f"bcftools norm -f {ref} -o {out_ln} {input}.gz")
        run_cmd(f"cat {out_ln} | sed 's/SOMATIC;//' | {vcf2table} > {out}")
        output_files.append(out)
        # cleanup
        shell(f"rm {input}.gz; rm {input}.gz.tbi; rm {input}; mv {out_ln} {input}")

else:
    for input in input_files:
        # print('input:', input)
        out = f"{input}.table"
        output_files.append(out)
        # varscan output has to be converted to avinput file format
        # anno_format adjusts alt, ref and start and end positions for indels
        run_cmd(f"{varscan2table} < {input} > {out}")

# CONCAT THE FILES
# rm first two lines (the headers of both files)
run_cmd(f"cat {' '.join(output_files)} | sort -V -k1,2 | mawk 'NR > 2 {{ print }}' >  {output}")

shell(f"rm -f {' '.join(output_files)}")
show_output(f"Concated {input_files[0]} and {input_files[1]} into {output} for annotation.", color="success")
# shell(f"rm {' '.join(output_files)}")
