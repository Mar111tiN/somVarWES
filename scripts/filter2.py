import os
from filters.filter2_core import get_filter2


def main(s):
    '''
    snakemake wrapper for calling the filter functions
    '''

    fconfig = s.config['filter']
    filter_file = os.path.join(
        s.config['paths']['filter_settings'],
        fconfig['filter_settings']
    )

    get_filter2(
        mut_file=s.input.filter1,
        filter2_output=s.output.filter2,
        filter_file=filter_file,
        filter_sheet=fconfig['excel_sheet'],
        filter_name=fconfig['filter2'],
        keep_syn=fconfig['keep_syn'],
        filterbam_output=s.output.filter2_for_filterbam,
        filterbam_stringency=s.config['filter_bam']['stringency_for_bam'],
        excel_output=fconfig['excel_output']
    )


if __name__ == "__main__":
    main(snakemake)
