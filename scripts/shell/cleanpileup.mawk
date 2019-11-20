#!/bin/sh

mawk '
{
    read = $0;
    # remove position traces
    gsub(/\^.|\$/,"",read);
    while (match(read,/[+-][0-9]+/)) {
        pre = substr(read,1,RSTART-2);
        indel = substr(read,RSTART,1); # + or -
        base = substr(read,RSTART-1,1); # base before the deletion
        l = substr(read,RSTART+1,RLENGTH-1)
        post = substr(read,RSTART+RLENGTH+l)
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
    print read;
}'