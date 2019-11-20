#!/bin/sh
mawk '
BEGIN {
    OFS="\t";
}
NR > 1 {
    start = $2 - 1;
    if ($5 == "-") {
        start = start -1;
    }
    print($1,start,start+1);
}'
