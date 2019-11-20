#!/bin/sh

mawk '
BEGIN {
    # ROWS --> DATA-fields
    Row["NR1"]="5";
    Row["NR2"]="6";
    Row["TR1"]="9";
    Row["TR2"]="10";
    Row["somatic_status"]="13";
    Row["NR1+"]="20";
    Row["NR2+"]="22";
    Row["TR1+"]="16";
    Row["TR2+"]="18";
    Row["somaticP"]="15";
    Row["variantP"]="14";

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
    printf("Chr\tStart\tEnd\tRef\tAlt");
    for (i = 4; i++ < 16;) {
        printf("\t%s",HEADER[i])
    }
    printf("\n");
}
NR > 1 {
    start=$2;
    end=$2; 
    # variant is an insertion
    if (sub(/\+/,"",$4)) {
        $3="-";
    # variant is a deletion
    } else if (sub(/\-/,"",$4)) {
        $3=$4;
        $4="-";
        start=$2+1;
        end=$2+length($3);
    }
# print all fields
    printf("%s\t%s\t%s\t%s\t%s",$1,start,end,$3,$4);
    for (i = 4; i++ < 15;) {
        printf("\t%s", $Row[HEADER[i]])
    }
    printf("\n");
}'