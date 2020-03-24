from os import system as shell
from script_utils import show_output, show_command

##############################################################################
# ######################### SNAKE PARAMETERS ##################################


def main(s):
    w = s.wildcards
    threads = s.threads
    # log = s.log
    matrix_file = s.output[0]
    params = s.params
    pon_list = params.pon_list
    EBparams = s.config['EBFilter']['params']
    bed_file = s.input[0]
    chrom = w.chrom
    i = w.i

    # import the scripts
    cleanpileup = params.cleanpileup
    pon2cols = params.pon2cols
    pile2count = params.pile2count

    # ############## PILEUP --> MATRIX FILE ##################
    # do the pileup into the matrix file

    show_output(f"Performing pileup and read matrix generation for split {i} of {chrom}", color='normal', time=True)
    pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']} -l {bed_file} -r {chrom} -b {pon_list}"
    # cut -f $({pon2cols}< {sample_list}) creates a cut command only including the desired
    pipe_cmd = f"{pileup_cmd} | cut -f $({pon2cols} < {pon_list}) | {cleanpileup} | {pile2count} | pigz -p {threads} > {matrix_file}"
    # do the pileup to matrix_file
    show_command(pipe_cmd, multi=False)
    shell(pipe_cmd)
    # cleanup if matrix file was generated
    shell(f"rm {bed_file}")
    show_output(
        f"Matrix generation for split {i} of {chrom} done! File saved to {matrix_file}",
        color='success',
        time=True
    )


if __name__ == "__main__":
    main(snakemake)
