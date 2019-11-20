#!/bin/sh
mawk '
BEGIN {
    OFS="\t
}
    {key=$1 "-" $2; print(key,$3,$4,$5}

}'
