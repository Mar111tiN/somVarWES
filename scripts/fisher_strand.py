from fisher_core import fisher


def main(s):
    fisher(s.input, s.output, threads=s.threads, log=s.log)


if __name__ == "__main__":
    main(snakemake)
