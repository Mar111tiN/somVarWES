from primer3_core import primer3_main


def main(s):
    # ############ SNAKEMAKE ##################

    config = s.config
    pconfig = config['primer3']

    PCR_config = {
        'seq_len': 500,
        'mut_pad': 5,
        'prod_size_min': pconfig['min_size'],
        'prod_size_max': pconfig['max_size']
    }
    primer3_main(
        s.input[0],
        s.output[0],
        PCR_config=PCR_config,
        chroms_folder=s.params.genome_split,
        threads=pconfig['threads']
        )


if __name__ == "__main__":
    main(snakemake)
