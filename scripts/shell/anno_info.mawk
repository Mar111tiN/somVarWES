#!/bin/sh

mawk '
NR == 1 {
    header = $0
    # DATA-fields --> HEADERS
    HEADER[5] = "somatic_status";
    HEADER[6] = "TR1";
    HEADER[7] = "TR1+";
    HEADER[8] = "TR2";
    HEADER[9] = "TR2+";
    HEADER[10] = "NR1";
    HEADER[11] = "NR1+";
    HEADER[12] = "NR2";
    HEADER[13] = "NR2+";
    HEADER[14] = "somaticP";
    HEADER[15] = "variantP";

    ######## HEADER #############
    expanded_info = HEADER[5]
    for (i = 5; i++ < 15;) {
        expanded_info = expanded_info "\t" HEADER[i]
    }
    sub("Otherinfo", expanded_info, header);
    print header;
}
NR > 1 {
    print $0
}'