from os import system as shell


def main(s):
    w = s.wildcards
    threads = s.threads
    config = s.config
    # log = s.log
    output = s.output

    if config['setup']['rerun']:
        rerun_bam = s.params.rerun_bam
        merge_cmd = f"ln -s {rerun_bam} {output.bam}; sambamba index {output.bam}"
        print(merge_cmd)
        shell(merge_cmd)
    else:
        # merge all split bams into one bam using sambamba merge (also sorts the merged bam)
        merge_cmd = f"sambamba merge -t {threads} {output.bam} bamfinalsplit/{w.sample}_{w.type}.*.bam; sambamba index {output.bam}"
        print(merge_cmd)
        shell(merge_cmd)


main(snakemake)