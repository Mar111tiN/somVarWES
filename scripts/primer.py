from primer3_core import primer3_main


def main(s):
    # ############ SNAKEMAKE ##################

    c = s.config
    cc = config["primer3"]

    PCR_config = {
        "seq_len": 500,
        "mut_pad": 5,
        "prod_size_min": cc["min_size"],
        "prod_size_max": cc["max_size"],
    }
    primer3_main(
        s.input[0],
        s.output[0],
        PCR_config=PCR_config,
        chroms_folder=s.params.genome_split,
        threads=cc["threads"],
    )


if __name__ == "__main__":
    main(snakemake)
