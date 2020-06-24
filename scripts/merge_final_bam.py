from os import system as shell

def main(s):
    w = s.wildcards
    threads = s.threads
    config = s.config
    # log = s.log
    output = s.output

    if config['setup']['rerun']:
        rerun_bam = s.params.rerun_bam
        shell(f"ln -s {rerun_bam} {output.bam}; sambamba index {output.bam}")
    else:
        # merge all split bams into one bam using samtools merge (also sorts the merged bam)
        shell(f"sambamba merge -t {threads} bamfinalsplit/{w.sample}_{w.type}.*.bam; ")


if __name__ == "__main__":
    main(snakemake)