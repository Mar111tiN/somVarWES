#!/bin/sh

mawk '
NR == 1 {
    header = $0
    # DATA-fields --> HEADERS
    HEADER[1] = "somatic_status";
    HEADER[2] = "TR1";
    HEADER[3] = "TR1+";
    HEADER[4] = "TR2";
    HEADER[5] = "TR2+";
    HEADER[6] = "NR1";
    HEADER[7] = "NR1+";
    HEADER[8] = "NR2";
    HEADER[9] = "NR2+";
    HEADER[10] = "somaticP";
    HEADER[11] = "variantP";

    ######## HEADER #############
    expanded_info = HEADER[5]
    for (i = 0; i++ < 11;) {
        sub("Otherinfo" i, HEADER[i], header);
    }
    
    print header;
}

NR > 1 {
    print $0
}'