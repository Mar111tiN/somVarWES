import os
import pandas as pd

# add ebscore package to sys.path
sys.path.append(os.path.join(snakemake.scriptdir, "myCNV/code"))
from pon_coverage import make_PON_coverage
from script_utils import show_output


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config
    p = s.params
    o = s.output
    input_list = s.input
    full_output_list = o.fullCov
    filter_output_list = o.filterCov

    # the sample list is that huge blob of files
    # samples should be processed chromosome-wise so chroms and samples have to be filtered out from input_list
    # easiest to be done using pandas :-)
    sample_df = pd.DataFrame(input_list, columns=["file"])
    sample_df[["sample", "Chr"]] = sample_df["file"].str.extract(
        r"([^/]+)\.(chr[0-9XY]+)\."
    )

    # load all sample coverages for one chromosome
    show_output(
        f"Loading all PON coverages from {p.pon_path} for exome-wide normalization",
        time=True,
    )
    # output the file
    show_output(f"Combining PON coverage for {len(sample_df.index)} samples", time=True)

    show_output(f"{len(sample_df.index)} samples detected:")
    full_df, filter_df = make_PON_coverage(sample_df, config=c["CNV"]["PONcoverage"])

    # take the output_file_name from first element of output_base
    full_output_base = full_output_list[0]
    filter_output_base = filter_output_list[0]

    for chrom in sample_df["Chr"].unique():
        full_file = full_output_base.replace("chr1", chrom)
        filter_file = filter_output_base.replace("chr1", chrom)

        show_output(
            f"Writing full coverage for chrom {chrom} to {full_file}", color="success"
        )
        full_df.query("Chr == @chrom").to_csv(
            full_file, sep="\t", index=False, compression="gzip"
        )

        show_output(
            f"Writing filtered coverage for chrom {chrom} to {filter_file}",
            color="success",
        )
        filter_df.query("Chr == @chrom").to_csv(
            filter_file, sep="\t", index=False, compression="gzip"
        )


if __name__ == "__main__":
    main(snakemake)
