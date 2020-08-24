#!/bin/sh

# specialized converter for HDRtools
# use samtools mpileup -f ! bam | pile2hotspot [MinALT=3]
# takes pileup data from stream and outputs mutations with
# ..with a minimum occurrence of MinAlt
 
mawk '
# get the MAPQ as argument with default 20
BEGIN {
    MinAlt = "'${1:-3}'" + 0;
    hCols="Chr,Pos,Ref,Depth,Alt,AltSum,AltRatio,totalAlt";
    hColsCount = split(hCols,HCOLS,",");
    for (col = 0; col++ < hColsCount-1;) {
        printf("%s\t",HCOLS[col]);
    }
    printf("%s\n",HCOLS[hColsCount]);

    ##### PREPARE LETTER COUNT #######
    # get the letters to look for
    letterCount = split("Aa,Gg,Cc,Tt,+,-",LETTERS,",")
    # LETTER
    ###### HEADER ################
    # create the patterns to look for
    for (l=0;l++<letterCount;) {
        LETPAT[l] = "[" LETTERS[l] "]"
    }
}
$5 ~ /[AaCcGgTt]/ {
    totalAlt=0;
    maxCount=0;
    maxLetter=0;
    # loop through the letters and count the occurrences
    for (l = 0;l++< letterCount;) {
        COUNT[l] = gsub(LETPAT[l], "", $5);
        totalAlt = totalAlt + COUNT[l];
        if (COUNT[l] > maxcount) {
            maxCount = COUNT[l];
            maxLetter = substr(LETTERS[l],1,1);
        }
    }
    # OUTPUT only above a minimum of MinAlt mutations #########
    if (totalAlt >= MinAlt) {
        # print base columns
        for (l=0;l++<4;) {
            printf("%s\t", $l);
        }
        altRatio=maxCount/$4;
        printf("%s\t%s\t%s\t%s\n",maxLetter,maxCount,altRatio,totalAlt);
    }
}'