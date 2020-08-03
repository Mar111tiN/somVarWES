from anno_core import run_annovar


def main(s):
    p = s.params
    run_annovar(s.input, s.output,
                annovar=s.config['tools']['annovar'],
                annoinfo=p.anno_info,
                annoparams=p.anno
                )


if __name__ == "__main__":
    main(snakemake)
