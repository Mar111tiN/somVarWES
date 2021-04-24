import os
from filters.filter2_core import get_filter2


def main(s):
    '''
    snakemake wrapper for calling the filter functions
    '''

    p = s.params
    cc = s.config["filter"]

    filter_config = dict(
        filter_file=os.path.join(p.snakedir, cc["filter_settings"]),
        filter_sheet=cc["excel_sheet"],
        filter_name=cc["filter1"],
        filterbam_stringency=s.config['filter_bam']['stringency_for_bam'],
        excel_output=cc['excel_output'],
        keep_syn=cc['keep_syn'],
        pop_cols=cc["pop_cols"]
    )

    get_filter2(
        mut_file=s.input.filter1,
        filter2_output=s.output.filter2,
        config=filter_config,
        filterbam_output=s.output.filter2_for_filterbam,
    )


if __name__ == "__main__":
    main(snakemake)
