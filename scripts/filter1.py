import os
from filters.filter1_core import filter_basic_filter1


def main(s):
    """
    snakemake wrapper for calling the filter functions
    """

    p = s.params

    cc = s.config["filter"]

    filter_config = dict(
        filter_file=os.path.join(p.snakedir, cc["filter_settings"]),
        filter_sheet=cc["excel_sheet"],
        filter_name=cc["filter1"],
        keep_syn=cc["keep_syn"],
        keep_UTR=cc["keep_UTR"],
        pop_cols=cc["pop_cols"],
    )

    filter_basic_filter1(
        mut_file=s.input[0],
        basic_output=s.output.basic,
        filter1_output=s.output.filter1,
        config=filter_config,
    )


if __name__ == "__main__":
    main(snakemake)
