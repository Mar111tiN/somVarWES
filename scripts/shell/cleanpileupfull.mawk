#!/bin/sh

### cleans samtools mpileup output 
# removes from all reads the traces of base location and indel lengths

mawk '
NR == 1 {
    samples = (NF - 3) / 3;
}
{
    read = $0;
    
    # remove position traces from all read fields
    for (i = 0; i++ < samples;) {
        col = (i*3) + 2;
        gsub(/\^[^\t]|\$/,"",$col);
    }
    # 
    while (match(read,/[+-][0-9]+/)) {
        pre = substr(read,1,RSTART-2);
        indel = substr(read,RSTART,1); # + or -
        base = substr(read,RSTART-1,1); # base before the deletion
        l = substr(read,RSTART+1,RLENGTH-1);
        post = substr(read,RSTART+RLENGTH+l);
        if (indel == "-") {
            if (match(base,/[ACGT]/)) {
                base = "D";
            } else {
                base = "d";
            }
        } else {
            if (match(base,/[ACGT]/)) {
                base = "I";
            } else {
                base = "i";
            }            
        }

        read = pre base post;
    }     
# print all fields
    for (i=0; i++ < NF;) {
        printf("%s\t",$i);
    }
    printf("\n");
}'
