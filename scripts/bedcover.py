from coverage import get_cover_svg

import os
import sys

def main(s):

    workdir = s.config['workdir']
    input = s.input
    output = os.path.join(workdir, str(s.output))

    params = s.params
    w = s.wildcards
    sample_name = f"{w.sample}{w.type}"

    get_cover_svg(input.sample, output, sample=sample_name,
                  exon_cover=params.exon_cover, refgen=params.refgen, log=s.log, prettifyBed=params.prettifyBed)


if __name__ == "__main__":
    main(snakemake)
