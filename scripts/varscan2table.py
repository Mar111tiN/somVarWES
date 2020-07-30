import os


def main(s):
    config = s.config

    params = s.params
    varscan2table = params.varscan2table

    build = config['ref']['build']
    ref = os.path.join(config['paths']['mystatic'],
                       config['ref'][build]['genome'])

    varscan2table(s.input, s.output,
                  refgen=ref,
                  isVCF=config['varscan']['vcf'],
                  vcf2csv=params.vcf2csv,
                  editcsv=params.editcsv,
                  coords2annovar=params.coords2annovar
                  )


if __name__ == "__main__":
    main(snakemake)
