import os
from HDR_core import run_HDR


def main(s):
    # ############## s ################################
    w = s.wildcards
    config = s.config
    i = s.input
    p = s.params
    PAD = min(config['HDR']['padding'], config['filter_bam']['padding'])
    bam = os.path.split(s.input.filter_bam)[1]
    bam = os.path.join(config['filter_bam']['folder'], bam)

    # get the tumor and normal bam files from the tumor-normal file name
    # remove the normal part of the tumor-normal descriptor
    tumor_bam = bam.replace(
        f"-{w.normal}", '').replace('.done', '.bam').replace('filterbamdone', 'filterbam')
    # remove the tumor part of the tumor-normal descriptor
    normal_bam = bam.replace(
        f"{w.tumor}-", '').replace('.done', '.bam').replace('filterbamdone', 'filterbam')
    print('tumor:', tumor_bam, 'normal:', normal_bam)

    # ## run the main HDR function
    run_HDR(
        mut_file=i.filter_file,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        filter_pileup=i.pileup,
        out_file=s.output.HDR_table,
        MINSIM=p.min_sim,
        PAD=PAD,
        MINQ=p.min_q,
        HDRMINCOUNT=p.min_HDR_count
    )


if __name__ == "__main__":
    main(snakemake)