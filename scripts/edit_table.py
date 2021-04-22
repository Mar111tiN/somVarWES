import os
import pandas as pd
from script_utils import sort_df
from derive_cols import get_PON_info, apply_clinscore, get_gene_lists, rearrange_cols


def main(s):
    """
    run the editing of the mutation_list
    """

    input = str(s.input)
    output = str(s.output)
    p = s.params
    config = s.config
    cc = config["editList"]
    snake_folder = p.snake_folder

    # ## LOAD MUT FILE #####################################
    print(f"Started editing and basic filtering for {input}.")
    anno_df = pd.read_csv(input, sep="\t")

    # ########### DERIVE PON COLUMNS ##########################
    anno_df = get_PON_info(anno_df)

    # ########### APPLY CLINSCORE ##########################
    clinscore_file = cc["clinscore_yaml"]
    # get the clinscore dicts from yaml file
    if not clinscore_file.startswith("/"):
        clinscore_file = os.path.join(snake_folder, clinscore_file)

    anno_df = apply_clinscore(anno_df, clinscore_file)

    # ########### APPLY GENE LISTS ##########################
    candidate_excel = os.path.join(snake_folder, cc["candidate_list"])
    anno_df = get_gene_lists(anno_df, candidate_excel=candidate_excel)

    # ############## REARRANGE COLUMNS  and sort rows ###################################
    info_folder = os.path.join(snake_folder, "info")
    anno_df = rearrange_cols(anno_df, file_name="filter_cols", folder=info_folder)
    anno_df = sort_df(anno_df)

    anno_df.to_csv(output, sep="\t", index=False)
    print(f"Writing edited mutation list with added columns to {output}.")


if __name__ == "__main__":
    main(snakemake)
