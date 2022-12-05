import os
import pandas as pd
from script_utils import sort_df, show_output
from derive_cols import (
    get_PON_info,
    apply_clinscore,
    get_gene_lists,
    rearrange_cols,
    get_refgene,
    addGenmap,
    addGCratio,
)


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
    genmap_folder = p.genmap_folder
    gc_folder = p.gc_folder

    cc["keep_syn"] = config["filter"]["keep_syn"]
    cc["keep_UTR"] = config["filter"]["keep_UTR"]

    # ## LOAD MUT FILE #####################################
    show_output(f"Started editing and basic filtering for {input}.")
    anno_df = pd.read_csv(input, sep="\t")

    # ########### change refGene COLUMNS ######################
    anno_df = get_refgene(anno_df, config=cc)
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
    # ########### ADD MAPPABILITY + GCratio ###########################
    anno_df = addGenmap(anno_df, genmap_path=genmap_folder, chrom_list=p.chrom_list)
    anno_df = addGCratio(anno_df, gc_path=gc_folder, chrom_list=p.chrom_list)
    # ############## REARRANGE COLUMNS  and sort rows ###################################
    info_folder = os.path.join(snake_folder, "info")
    anno_df = rearrange_cols(anno_df, file_name="filter_cols", folder=info_folder)
    anno_df = sort_df(anno_df)
    anno_df.loc[:, "End"] = anno_df["End"].astype(int)

    anno_df.drop_duplicates().to_csv(output, sep="\t", index=False)
    show_output(
        f"Writing edited mutation list with added columns to {output}.", color="success"
    )


if __name__ == "__main__":
    main(snakemake)
