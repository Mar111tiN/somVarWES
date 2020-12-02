#!/bin/sh

mawk '
# KEEP HEADER
NR == 1 {
    print $0;
    # get samplepos as bash arg
    samplepos='"$1"';
}
NR > 1 {
    # print constant fields
    printf("%s", $1)
    for (i = 1; i++ <5;) {
       printf("\t%s", $i); 
    }
    # go through fields
    for (i = 5; i++ < NF;) {
        # split data by "|"
        split($i,array,"|");
        count=0;
        for (j = 0; j++ < length(array);) {
            # get the position argument
            if (j != samplepos) {
                count++;
                newarray[count]=array[j];
            } else {
                target=array[j];
            }
        }
        #output
        printf("\t%s", target);
        for (r=0; r++ < length(newarray);) {
            printf("|%s", newarray[r]);
        }
    }
    printf("\n");
}
'