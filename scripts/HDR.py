import os
from HDR_core import HDR_master


def main(s):
    # ############## s ################################
    w = s.wildcards
    config = s.config
    i = s.input
    p = s.params
    threads = config['HDR']['threads']
    bam = i.bam

    # load all HDR params into HDR_config
    #
    HDR_config = {
        'MINSIM': config['HDR']['min_similarity'],
        'MINQ': config['mpileup']['Q'],
        'MINq': config['mpileup']['MAPQ'],
        'MinHDRCount': config['HDR']['min_HDR_count'],
        'PAD': min(config['HDR']['padding'], config['filter_bam']['padding']),
        'MinAltSupport': config['HDR']['min_alt_support']
    }

    # add the editbamdf tool to the config
    HDR_config['editbamdf'] = p.editbamdf
    # ## run the main HDR function

    HDR_df = HDR_master(
        mut_file=i.filter_file,
        bam_file=i.bam,
        chrom=w.chrom,
        pileup_file=i.pileup,
        threads=threads,
        HDR_config=HDR_config,
        out_file=s.output.HDR_table
    )


if __name__ == "__main__":
    main(snakemake)
