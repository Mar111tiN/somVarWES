import os
from HDR_run import HDR_master
from script_utils import show_output


def main(s):
    # ############## s ################################
    w = s.wildcards
    config = s.config
    i = s.input
    p = s.params
    threads = config['HDR']['threads']
    bam_file = i.bam
    out_file = s.output.HDR_table
    # load all HDR params into HDR_config
    #
    HDR_config = {
        "minAltSum": config['HDR']['min_alt_sum'],
        "minAltRatio": config['HDR']['min_alt_ratio'],
        "maxAltRatio": config['HDR']['max_alt_ratio'],
        'MINSIM': config['HDR']['min_similarity'],
        'MINQ': config['mpileup']['Q'],
        'MINq': config['mpileup']['MAPQ'],
        'MinHDRCount': config['HDR']['min_HDR_count'],
        'PAD': min(config['HDR']['padding'], config['filter_bam']['padding']),
        'MinAltSupport': config['HDR']['min_alt_support'],

    }

    # add the genomes folder to config
    HDR_config['genome_split_path']: p.genome_split
    # add the tools to the config
    HDR_config['editbam'] = p.editbam
    HDR_config['bam2csv'] = p.bam2csv
    HDR_config['pile2hotspot'] = p.pile2hotspot
    HDR_config['pile2hotspot_chrom'] = p.pile2hotspot_chrom
    # ## run the main HDR function

    HDR_df = HDR_master(
        mut_file=i.filter_file,
        bam_file=bam_file,
        chrom=w.chrom,
        threads=threads,
        HDR_config=HDR_config,
        pileup_file=i.pileup
    )

    HDR_df.to_csv(out_file, sep='\t', index=False)
    show_output(
        f"HDRdetect of {bam_file} finished. Writing results to {out_file}", time=True, color='success'
    )


if __name__ == "__main__":
    main(snakemake)
