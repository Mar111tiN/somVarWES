#!/bin/sh

mawk '
END {
    printf("%s,%s",1,2);
    for (i=0; i++<NR;) {
        printf(",%s",2+3*i)
    }
}'