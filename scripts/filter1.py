import os
from filters.filter1_core import filter1, filter1_AML7


def main(s):
    '''
    snakemake wrapper for calling the filter functions
    '''

    fconfig = s.config['filter']
    filter_file = os.path.join(
        s.config['paths']['filter_settings'],
        fconfig['filter_settings']
    )

    # use the right filter function
    if fconfig['filter1'] == "filter1":
        filter1(
            mut_file=s.input[0],
            basic_output=s.output.basic,
            filter1_output=s.output.filter1,
            filter_file=filter_file,
            filter_sheet=fconfig['excel_sheet'],
            filter_name=fconfig['filter1'],
            keep_syn=fconfig['keep_syn']
        )

    elif fconfig['filter1'] == "filter1-AML7":
        filter1_AML7(
            mut_file=s.input[0],
            basic_output=s.output.basic,
            filter1_output=s.output.filter1,
            filter_file=filter_file,
            filter_sheet=fconfig['excel_sheet'],
            filter_name=fconfig['filter1'],
            keep_syn=fconfig['keep_syn']
        )


if __name__ == "__main__":
    main(snakemake)
