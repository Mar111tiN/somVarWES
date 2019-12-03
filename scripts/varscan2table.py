from os import system as shell
import os
w = snakemake.wildcards
config = snakemake.config


input_files = snakemake.input
output = snakemake.output
params = snakemake.params
VCF = config['varscan']['vcf']
vcf2table = params.vcf2table
varscan2table = params.varscan2table
log = snakemake.log
build = config['ref']['build']
ref = os.path.join(config['paths']['mystatic'], config['ref'][build]['genome'])


# anno_format = params.anno_format
# expand_info = params.expand_info
# merge_anno = params.merge_anno
output_files = []
if VCF:  
    for input in input_files:
        # standard output
        print('input:', input)
        out = input.replace('vcf', 'table')
        out_ln = input.replace('vcf', 'ln.vcf')
        # for indel realignment, files have to be bgzipped and indexed with tabix
        shell(f"bgzip < {input} > {input}.gz")
        shell(f"tabix {input}.gz")
        shell(f"bcftools norm -f {ref} -o {out_ln} {input}.gz")
        shell(f"cat {out_ln} | sed 's/SOMATIC;//' | {vcf2table} > {out}")
        output_files.append(out)
        # cleanup
        shell(f"rm {input}.gz; rm {input}.gz.tbi; rm {input}; mv {out_ln} {input}")

else:
    for input in input_files:
        print('input:', input)
        out = f"{input}.table"
        output_files.append(out)
        # varscan output has to be converted to avinput file format
        # anno_format adjusts alt, ref and start and end positions for indels
        shell(f"{varscan2table} < {input} > {out}")

# CONCAT THE FILES

########### !! does not work for empty files !! ### --> workaround:
# no header and merge files with fixed varscan header in merge_table


# # get first line
# head_cmd = f"mawk 'NR == 1 {{print}}' < {output_files[0]} > {output}"
# print(head_cmd)
# shell(head_cmd)
# # concat and sort the files and append to output
# cat_cmd = f"cat {' '.join(output_files)} | sort -V -k1,2 | mawk 'NR > 2 {{ print }}' >>  {output}"
# print(cat_cmd)
# shell(cat_cmd)
####################################################################

# concat and sort the files and write to output
cat_cmd = f"cat {' '.join(output_files)} | sort -V -k1,2 | mawk 'NR > 2 {{ print }}' >  {output}"
print(cat_cmd)
shell(cat_cmd)

shell(f"rm -f {' '.join(output_files)}")
print(f"Concated {input_files[0]} and {input_files[1]} into {output} for annotation.")
# shell(f"rm {' '.join(output_files)}")
