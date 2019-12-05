import os
from os import system as shell
from datetime import datetime as dt


###########################################################
###################### FUNCTIONS ##########################


ansii_colors = {
          'magenta': '[1;35;2m',
          'green': '[1;9;2m',
          'red': '[1;31;1m',
          'cyan': '[1;36;1m',
          'gray': '[1;30;1m',
          'black': '[0m'
          }

colors = {
        'process': ansii_colors['green'],
        'time': ansii_colors['magenta'],
        'normal': ansii_colors['gray'],
        'warning': ansii_colors['red'],
        'success': ansii_colors['cyan']
        }

def show_output(text, color='normal', multi=False, time=False):
    '''
    get colored output to the terminal
    '''
    time = f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ''
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ''
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text)


def show_command(command, list=False, multi=True):
    '''
    prints the command line if debugging is active
    '''

    proc = f"\033[92mProcess {os.getpid()}\033[0m : " if multi else ""
    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(proc + command)
    return


##############################################################################
########################## SNAKE PARAMETERS ##################################

w = snakemake.wildcards
config = snakemake.config
threads = snakemake.threads
log = snakemake.log
matrix_file = snakemake.output[0]
params = snakemake.params
pon_list = params.pon_list
EBparams = snakemake.config['EBFilter']['params']
bed_file = snakemake.input[0]
chrom = w.chrom
i = w.i

# import the scripts
cleanpileup = params.cleanpileup
pon2cols = params.pon2cols
pile2count = params.pile2count

#############################################################################
############################### RUN #########################################

############### PILEUP --> MATRIX FILE ##################
base_file = bed_file.replace(".bed", "").replace("bed/", "")
#### !!!!! could also run the whole process multithreaded
#### have to use processor_id for files --> complicated

# do the pileup into the matrix file

show_output(f"Performing pileup and read matrix generation for split {i} of {chrom}", color='normal', time=True)
pileup_cmd = f"samtools mpileup -B -q {EBparams['MAPQ']} -Q {EBparams['Q']} -l {bed_file} -r {chrom} -b {pon_list}"
# cut -f $({pon2cols}< {sample_list}) creates a cut command only including the desired
pipe_cmd = f"{pileup_cmd} | cut -f $({pon2cols} < {pon_list}) | {cleanpileup} | {pile2count} | gzip > {matrix_file}"
# do the pileup to matrix_file
show_command(pipe_cmd, multi=False)
shell(pipe_cmd)
# cleanup if matrix file was generated
shell(f"rm {bed_file}")
show_output(f"Matrix generation for split {i} of {chrom} done! File saved to {matrix_file}", color='success', time=True)
