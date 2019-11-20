#!/bin/sh
mawk '
NR==1 { 
    poncount = NF - 2;

    ###### QUERY ##############
    # get the letters to look for
    split("AaGgCcTtIiDd",Letters,"")
    # get the length before depth addition
    len = length(Letters)
    Letters[len+1] = "Depth-ACGT";
    Letters[len+2] = "Depth-acgt";
    Letters[len+3] = "Depth-INDEL";
    Letters[len+4] = "Depth-indel";
    ###### HEADER ################

    printf("Chr\tStart");
    ####### HAVE TO ADJUST ###########
    # loop through the letters
    for (i = 0; i++ < len+4;) {
        printf("\t%s", Letters[i]);
    }
    printf("\n");
}

{   ######### LINES #############
    printf("%s\t%s",$1,$2);
   
    # loop through the reads
    for (r = 2; r++ < poncount + 2;) {
        read = $r;
        ######### --> COUNTS ############
        # loop through the letters
        for (i = 0; i++ < len;) {
            Count[r "-" i] = gsub(Letters[i], "", read);
        }
        ######## --> DEPTHS ##############
        # reset the depths
        for (i = 0; i++ < 4;){
            Count[r "-" len+i] = 0
        }  
        # Depth-ACGT - 1|3|5|7
        for (i = 1; i < 9; i+=2) {         
            Count[r "-" len+1] += Count[r "-" i];
        }
        # Depth-acgt - 2|4|6|8
        for (i = 2; i < 9; i+=2) {
            Count[r "-" len+2] += Count[r "-" i];
        }
        # Depth-INDEL
        for (i = 9; i < 12; i+=2) {
            Count[r "-" len+3] += Count[r "-" i];
        }
        # Depth-INDEL
        for (i = 10; i < 13; i+=2) {
            Count[r "-" len+4] += Count[r "-" i];
        }

    }
    ######### OUTPUT #############
    # loop through the letters
    printf("\t")
    for (i = 0; i++ < length(Letters);) {
        # first line extra for pretty
        printf("%s", Count[3 "-" i])
        for (r = 3; r++ < poncount + 2;) {
           printf("|%s", Count[r "-" i]);
        }
        printf("\t")  
    }
    printf("\n");
}'
