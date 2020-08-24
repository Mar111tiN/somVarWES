#!/bin/sh

# specialized converter for HDRtools
# use samtools mpileup -f ! bam | pile2hotspot chrom [MinALT=3]
# takes pileup data from stream and outputs mutations with
# ..a minimum occurrence of MinAlt
# derived from pile2hotspot but used for one chrom
# stops after given chrom has been read through

mawk '
# get the MAPQ as argument with default 20
BEGIN {
    chrom = "'$1'"
    if (chrom !~ "chr") {
        print "Please provide a chromosome as arg1!" > "/dev/stderr";
        exit;
    }
    MinAlt = "'${2:-3}'" + 0;
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
$1 !~ chrom {
    if (found == 1) exit;
    next;
    }
$5 ~ /[AaCcGgTt]/ {
    found=1;
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