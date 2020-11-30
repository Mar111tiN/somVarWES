from anno_core import run_annovar
import os

#################################
# add information to table/{sample}_{tumor_normal}.csv
def get_anno_params(c):
    """
    helper function to create full annovar parameters from cs
    """

    # shortcut to anno config
    cc = c["annovar"]
    # get the full path to humandb
    humandb = os.path.join(c["paths"]["mystatic"], c["annovar"]["humandb"])
    # get the available anno files
    file_list = list(os.walk(humandb))[0][-1]
    # reduce anno files to the files compatible with genome build version
    build = c["ref"]["build"]
    build_files = []
    for file in file_list:
        if build in file:
            build_files.append(file)

    # filter the anno protocol for the available files for that genome build
    anno_refs = cc["annofiles"]
    anno_list = []
    missing_list = []
    for anno in anno_refs:
        for file in build_files:
            if anno in file:
                anno_list.append(anno)
                break
        # if anno has not been found in file during for-loop
        else:
            missing_list.append(anno)

    # create the protocol string
    protocol = ",".join(anno_list)
    print(f"{' '.join(missing_list)} not found for {build}! Doing without.. ")
    # create the operation string 'g,r,f,f,f,f' assuming all but the first three dbs (ref, cytoBand, superDups) in config to be filter-based
    operation_list = []
    for anno in anno_list:
        if "Gene" in anno:
            operation_list.append("g")
        elif anno in ["cytoBand", "genomicSuperDups"]:
            operation_list.append("r")
        else:
            operation_list.append("f")
    operation = ",".join(operation_list)

    options = f'{humandb}/ -buildver {build} -remove -thread {cc["threads"]} -protocol {protocol} -operation {operation} -nastring "." -otherinfo'
    return options


def main(s):
    p = s.params
    run_annovar(
        s.input,
        s.output,
        annovar=s.config["tools"]["annovar"],
        annoinfo=p.anno_info,
        annoparams=get_anno_params(s.config),
    )


if __name__ == "__main__":
    main(snakemake)
