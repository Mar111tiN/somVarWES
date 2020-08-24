from io import StringIO
from subprocess import Popen, PIPE, run
import pandas as pd


def bam2df2(bam_file, HDR_config, bamtags="", tool='samtools', mut_row=None, region='', ):
    '''
    set the region requires 3 threads
    '''

    bam2csv = HDR_config['bam2csv']
    editbam = HDR_config['editbam']
    if region:
        mut_pos = region
    elif isinstance(mut_row, pd.Series):
        chrom = mut_row['Chr']
        pos = mut_row['Pos']
        mut_pos = f"{chrom}:{pos}-{pos} "
        print(mut_pos)
    else:
        mut_pos = ''
    cmd = f"{tool} view {bam_file} {mut_pos} | {bam2csv} | {editbam} {HDR_config['MINq']}"
    # show_command(cmd)
    bam_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return bam_df
