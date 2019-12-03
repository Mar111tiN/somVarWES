#!/bin/sh

mawk '
NR > 1 {
# snp_rwrt functionality
    start=$2;
    end=$2; 
# variant is an insertion
    if (sub(/\+/,"",$4)) {
        $3="-";
    }
# variant is a deletion
    if (sub(/\-/,"",$4)) {
        $3=$4;
        $4="-";
        start=$2+1;
        end=$2+length($3);
    }
# print all fields
    printf("%s\t%s\t%s\t",$1,start,end);
# format probabilities
    for (i=3; i<NF;i++) {
        if (match($i,"\.[0-9][0-9][0-9][0-9]")) {
            printf("%8.4e\t",$i);
        } else {
            printf("%s\t",$i);            
        }
    }   
    print $NF; 
}'
