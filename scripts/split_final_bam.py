from os import system as shell
import os

def main(s):
    w = s.wildcards
    threads = s.threads
    config = s.config
    # log = s.log
    o = s.output
    i = s.input

    if config['setup']['rerun']:
        # split the rerun bam directly into the bamfinalsplit folder
        rerun_bam = s.params.rerun_bam
        shell(f"sambamba view -t {threads} -h -f bam -o {o.bam} {rerun_bam} {w.chrom} ")
    else:
        # create symlinks from the proper chrom-split bams into bamfinalsplit
        print(f"ln -s $(pwd)/{i.bam} {o.bam}; sambamba index {o.bam}")
        shell(f"ln -s $(pwd)/{i.bam} {o.bam}; sambamba index {o.bam}")

main(snakemake)