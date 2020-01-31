import pandas as pd
import os
from subprocess import check_call as shell
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn
mpl.use('Agg')
seaborn.set()


def get_cover_svg(csv_file, sample):
    '''
    creates an svg-output file from the coverage summary
    '''

    data = pd.read_csv(csv_file, sep='\t', header=None)
    # sort out the necessary rows and remove first row (no coverage)
    cov = data.iloc[1:, [1, 2, 4, 5, 6, 7]]
    # add the respective columns
    cov.columns = ["coverage", "bases_at_coverage", "percentage_at_depth", "bases_at_min_coverage", "base_freq_at_min_coverage", "total_on_target"]
    # normalize depth percentage for simultaneous output into svg
    cov['percentage_at_depth'] = cov['percentage_at_depth'] / cov['percentage_at_depth'].max() * 100
    # calculate 95% coverage
    _95 = cov.query('base_freq_at_min_coverage > 95')['base_freq_at_min_coverage'].idxmin()

    plt.plot(cov['coverage'],cov['base_freq_at_min_coverage'], label='% of reads at coverage')
    plt.plot(cov['coverage'], cov['percentage_at_depth'], label='coverage distribution')
    plt.xlim(0, cov['percentage_at_depth'][10:].idxmax() * 3)
    plt.title(f'{sample} - 95% of reads over {_95}-fold coverage')
    plt.legend()
    plt.xlabel('coverage depth')
    # bla.text --> bla.svg
    svg_file = csv_file.replace('txt', 'svg')
    print(f"Saving {svg_file}")
    plt.savefig(f"{svg_file}")


workdir = snakemake.config['workdir']
input = snakemake.input
output = os.path.join(workdir, str(snakemake.output))
log = snakemake.log

params = snakemake.params
exon_cover = params.exon_cover
format_coverage = params.format_coverage
prettifyBed = params.prettifyBed
w = snakemake.wildcards
sample_name = f"{w.sample}{w.type}"

fastq_pair = ' '.join([os.path.join(workdir, fastq) for fastq in params.fastq])
refgen = params.refgen
histo_steps = snakemake.config['cover_bed']['histo_steps']
cmd = f"bedtools coverage -b {input.sample} -a {exon_cover} -hist -sorted -g {refgen}.genome 2>{log} | grep \'^all\' | sort -k2,2nr | {prettifyBed} | sort -k2,2n > {output}"
exit = shell(cmd, shell=True)
if exit == 0:
    get_cover_svg(output, sample_name)
# I don' think this is still needed - cover svg is much nicer
# cmd = f"{format_coverage} {output} {fastq_pair} {histo_steps} 1>{output}.summary 2>/dev/null"
shell(cmd, shell=True)
