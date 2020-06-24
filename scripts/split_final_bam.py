from os import system as shell


def main(s):
    w = s.wildcards
    threads = s.threads
    config = s.config
    # log = s.log
    output = s.output

    if config['setup']['rerun']:
        # split the rerun bam directly into the bamfinalsplit folder
        rerun_bam = s.params.rerun_bam
        shell(f"sambamba view -t {threads} -h -f bam -o {output.bam} {rerun_bam} {w.chrom} ")
    else:
        # create symlinks from the proper chrom-split bams into bamfinalsplit
        shell(f"ln -s {input.bam} {output.bam}; sambamba index {output.bam}")

main(snakemake)