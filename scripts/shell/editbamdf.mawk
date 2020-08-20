#!/bin/sh

# specialized bamcsv converter for HDRtools
# use samtools | bam2csv | editbam
 
mawk '
BEGIN {
    # get the MAPQ as argument with default 20 and convert to number
    MQ = "'${1:-20}'" + 0;
}
NR == 1 {
    # printf("MAPQ=%s\n", MQ);
    # define the output cols with their positions in the csv
    hCols="Chr,Pos,Seq,Qual,read_len,Soft_start";
    hColsCount = split(hCols,HCOLS,",");

    dCols="3,4,8,9"
    dColsCount = split(dCols,DCOLS,",");
    ######## HEADER OUTPUT #############
    for (col = 0; col++ < hColsCount-1;) {
        printf("%s\t",HCOLS[col]);
    }
    printf("%s\n",HCOLS[hColsCount]);
}

##### PROCESS LINES #################
NR > 1 && $5 >= MQ{
    read_len = length($8);
    if (match($6, "^[0-9]+S")) {
        soft_start = substr($6,RSTART,RLENGTH-1);
    } else {
        soft_start = 0;
    } 
    ######## OUTPUT #################
    for (col = 0; col++ < dColsCount;) {
        printf("%s\t", $DCOLS[col]);
    }
    printf("%s\t%s\n", read_len, soft_start)
}'