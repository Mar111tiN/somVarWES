import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn
from script_utils import run_cmd, show_output
mpl.use('Agg')
seaborn.set()


def make_svg(csv_file, plot_name="coverage_plot"):
    '''
    creates an svg-output file from the coverage summary
    '''

    data = pd.read_csv(csv_file, sep='\t', header=None)
    # sort out the necessary rows and remove first row (no coverage)
    cov = data.iloc[1:, [1, 2, 4, 5, 6, 7]]
    # add the respective columns
    cov.columns = ["coverage", "bases_at_coverage", "percentage_at_depth",
                   "bases_at_min_coverage", "base_freq_at_min_coverage", "total_on_target"]
    # normalize depth percentage for simultaneous output into svg
    cov['percentage_at_depth'] = cov['percentage_at_depth'] / \
        cov['percentage_at_depth'].max() * 100
    # calculate 95% coverage
    over95 = cov.query('base_freq_at_min_coverage > 0.95')
    if len(over95):
        _95 = over95['base_freq_at_min_coverage'].idxmin()
    else:
        _95 ="n.a"

    plt.plot(cov['coverage'], cov['base_freq_at_min_coverage'],
             label='% of reads at coverage')
    plt.plot(cov['coverage'], cov['percentage_at_depth'],
             label='coverage distribution')
    plt.xlim(0, cov['percentage_at_depth'][10:].idxmax() * 3)
    plt.title(f'{plot_name} - 95% of reads over {_95}-fold coverage')
    plt.legend()
    plt.xlabel('coverage depth')
    # bla.text --> bla.svg
    svg_file = csv_file.replace('txt', 'svg')
    print(f"Saving {svg_file}")
    plt.savefig(f"{svg_file}")


def get_cover_svg(i, o, sample, exon_cover, refgen, log, prettifyBed):
    cmd = f"bedtools coverage -b {i} -a {exon_cover} -hist -sorted -g {refgen}.genome | grep \'^all\' | sort -k2,2nr | {prettifyBed} | sort -k2,2n > {o}"
    if run_cmd(cmd):
        make_svg(o, plot_name=sample)
    else:
        show_output(f"bedtools for {i} exited with error!", color="warning")
