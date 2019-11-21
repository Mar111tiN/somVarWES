#!/bin/sh
# works like matrix2EBinput but treats inserts and deletions separately

input_fastq=$1;
splitFactor=$2;

gunzip < $1 | mawk '
BEGIN {
    file_base = "'$3'";
    splits = '$splitFactor'*4;
    for (x = 0; x++ < splits;) {
        toFile[x%splits]= int((x-1) / 4);
    }  
}
{
    file_number=toFile[NR%splits];
    print >> file_base"."file_number".fastq";
}

END {
    printf("NR=%s:file%s.\n", x, toFile[x]);
}'
