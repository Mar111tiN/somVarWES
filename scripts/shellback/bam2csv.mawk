#!/bin/sh

# transparent shell wrapper around a powerfull and fast bam translator
# shell layer allows argument passing
# maybe include help field
# only optional argument is a comma-separated list of additional tags 
# for 10x: CB,CY,UB,UY (default)
# for UMI-tools: RX,QX,MC,MD,RG,NM,AS,XS,MI,cD,cE,cd,ce
 
mawk '
BEGIN {
    ################ BAM2CSV TRANSLATOR ###################

    # define the header columns

    # here go the standard 11 tab-separated fields as described in the SAM-specs
    bamspex="QNAME,FLAG,CHR,START,MAPQ,CIGAR,RNEXT,LENGTH,PNEXT,SEQ,QUAL";
    # get the SPEX array to get the data from
    # SPEX[1] = "QNAME"
    # SPEX[2] = "FLAG"
    # --
    # SPEX[11] = "QUAL"
    spexCount = split(bamspex,SPEX,",");
    # here come the actual output fields. omitting RNEXT and PNEXT
    # can be flexibly changed
    outFields="QNAME,FLAG,CHR,START,MAPQ,CIGAR,LENGTH,SEQ,QUAL";
    ######### TAGS #####################
    # get extra tags from command argument $1 (with default for 10x)
    extra="'${1:-CB,CY,UB,UY}'";
    # and paste with the outFields
    cols = outFields "," extra

    # split the cols string into the DATAFIELDS array
    fieldCount = split(cols,FIELDS,",");
    # FIELDS[1] = "QNAME"
    # ..
    # FIELDS[9] = "QUAL"
    # FIELDS[10] = "CB"
    # ..
    # reverse the array for easy check
    for (i = 0; i++ < fieldCount;) {
      COLS[FIELDS[i]] = i;
    }
    # COLS["QNAME"] = 1 ....
    # COLS["CB"] = 10

    ######## HEADER OUTPUT #############
    for (col = 0; col++ < fieldCount-1;) {
      # now I can check for cols
      if (FIELDS[col] in COLS) { # is this neccessary
        printf("%s\t",FIELDS[col]);
      }
    }
    printf("%s\n",FIELDS[fieldCount]);
}

##### PROCESS LINES #################
{
    # get normal fields with COLS check
    for (col = 0; col++ < spexCount;) {
      if (SPEX[col] in COLS) {
        Data[SPEX[col]]=$col;
      }
    }
    # get TAGS 
    for (col = spexCount-2; col++ < fieldCount;) {
      # create the pattern from the DataField and the value format field  (like CB:Z:ATTGCCTG)
      pattern=FIELDS[col] ":[ZiB]:[^\\t]+\\t?[$]?";
      if (match($0, pattern)) {
        Data[FIELDS[col]]=substr($0,RSTART+5,RLENGTH-6);
      } else {
        Data[FIELDS[col]]=".";
      }
    }

    ######## OUTPUT #################
    for (col = 0; col++ < fieldCount - 1;) {
      if (FIELDS[col] in COLS) {
        printf("%s\t",Data[FIELDS[col]]);
      }
    }
    printf("%s\n",Data[FIELDS[fieldCount]]);
}'