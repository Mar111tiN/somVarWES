from os.path import getsize as filesize
from ebutils import matrix2AB_multi
from script_utils import show_output


def main(s):

    w = s.wildcards
    threads = s.threads
    # log = s.log
    cache_file = s.output[0]
    EBparams = s.config['EBFilter']['params']
    penalty = EBparams['fitting_penalty']
    matrix_file = s.input[0]
    chrom = w.chrom
    i = w.i

    # ############## MATRIX FILE --> EBcache ##################

    show_output(f"Performing AB-parameter computation for split {i} of {chrom}", color='normal', time=True)
    if filesize(matrix_file) > 20:
        cache_file = matrix2AB_multi(matrix_file, cache_file, penalty, threads)
        show_output(
            f"EBcache for split {i} of {chrom} done! Gzipped cache saved to {cache_file}",
            color='success',
            time=True
        )
    else:
        open(cache_file, 'a').close()
        show_output(
            f"Matrix file {matrix_file} was empty! Empty file {cache_file} is touched lest the pipeline break.",
            color='success',
            time=True
        )


if __name__ == "__main__":
    main(snakemake)
