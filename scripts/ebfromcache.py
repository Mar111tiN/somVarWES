import os
from ebrun import run_eb_from_cache


def main(s):
    w = s.wildcards
    p = s.params
    static_path = s.config['paths']['mystatic']
    EBconfig = s.config['EBFilter']
    pon_list = os.path.join(static_path, EBconfig['pon_list'])

    run_eb_from_cache(
        table=s.input.table,
        tumor_bam=s.input.tumor_bam,
        AB_cache_file=s.input.EBcache,
        matrix_cache_file=s.input.EBcache.replace('.cache', '.matrix'),
        output=s.output,
        pon_list=pon_list,
        chrom=w.chrom,
        log=s.log,
        threads=s.threads,
        EBparams=EBconfig['params'],
        full_output=EBconfig['full_pon_output'],
        # import the scripts
        cleanpileup=p.cleanpileup,
        csv2bed=p.csv2bed,
        pile2count=p.pile2count,
        matrix2EBinput=p.matrix2EBinput,
        reorder_matrix=p.reorder_matrix
    )


if __name__ == "__main__":
    main(snakemake)
