#!/fast/users/szyskam_c/work/miniconda/envs/hWES-env/bin/mawk -f

BEGIN {
# print the header
    printf("\n\n______________:::: Conda tools ::::_________________\n\n")
}
/perl/||/awk/||/bwa/||/gatk4/||/bedtools/||/bpipe/||/fastqc/||/picard/||/pindel/||/samtools/||/varscan/ {
    printf("  %-24s %-12s %12s\n", $1,$2,$4)
}
END {
    printf("\n___:::: paths / symlinks to specific versions ::::___\n \n")
}