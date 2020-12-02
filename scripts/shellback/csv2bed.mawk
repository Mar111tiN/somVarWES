#!/bin/sh
mawk '
BEGIN {
    OFS="\t";
}
$0 !~ "Chr" {
    start = $2 - 1;
    if ($5 == "-") {
        start = start -1;
    }
    print($1,start,start+1);
}'
