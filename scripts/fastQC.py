from os import system as shell

config = snakemake.config
w = snakemake.wildcards
path = f'{w.sample}_{w.type}_{w.readtrim}'
lines = round(config['qc']['samplefactor']) * 4
input = snakemake.input
log = snakemake.log

#   run QC on subsamples only if sample factor > 1
if config['qc']['samplefactor'] < 2:
    print(f"Performing fastQC for {input}")
    shell(f"fastqc {input} -o fastqc/ ") # &>{snakemake.log}
else:
    print(f"Downsampled for QC with factor: {config['qc']['samplefactor']}")
    shell(f"zcat {input} | mawk 'NR % {lines} > 0 && NR % {lines} < 5' > {input}.sub")
    shell(f'fastqc {input}.sub -o fastqc/ ') # &>{log}
    shell(f'mv fastqc/{path}.fastq.gz.sub_fastqc.zip fastqc/{path}_fastqc.zip')