#!/bin/sh

mawk '
BEGIN {
    # get the letters to look for
    split("AaGgCcTtIiDd",Letters,"")
    # get the length before depth addition
    len = length(Letters)
    Letters[len+1] = "Depth-ACGT";
    Letters[len+2] = "Depth-acgt";
    Letters[len+3] = "Depth-INDEL";
    Letters[len+4] = "Depth-indel";
    # get positions from letters
    # shift 5 columns for Chr...Alt
    for (x=0; x++ < len+4;) {
        letterPos[Letters[x]] = 5 + x;
    }
    printf("Chr\tStart\tEnd\tRef\tAlt\tdepthP\tmisP\tdepthN\tmisN\n")
}
NR > 1 {
    printf("%s\t%s\t%s\t%s\t%s\t", $1,$2,$3,$4,$5)
    Ref = $4;
    Alt = $5;
    if ((Ref == "-") || (Alt == "-")) {
        # get the values of the respective 
        split($letterPos["Depth-ACGT"],depthP,"|");
        split($letterPos["Depth-INDEL"],INDEL,"|");    
        split($letterPos["Depth-acgt"],depthN,"|");
        split($letterPos["Depth-indel"],indel,"|");
        split($letterPos["I"],I,"|");
        split($letterPos["D"],D,"|");    
        split($letterPos["i"],i,"|");
        split($letterPos["d"],d,"|");
        ############### FIRST VALUE ###############
        depth_p = depthP[1] + INDEL[1];
        mis_p = I[1] + D[1];
        depth_n = depthN[1] + indel[1];
        mis_n = i[1] + d[1];
        ############## LOOP THROUGH OTHER VALUES ####
        for (x = 1; x++ < length(depthP);) {
            depth_p = depth_p "|" (depthP[x] + INDEL[x]);
            mis_p = mis_p "|" (I[x] + D[x]);
            depth_n = depth_n "|" (depthN[x] + indel[x]);
            mis_n = mis_n "|" (i[x] + d[x]);
        }
    } else {
        depth_p = $letterPos["Depth-ACGT"];  
        mis_p = $letterPos[Alt];
        depth_n = $(letterPos["Depth-ACGT"]+1);  
        mis_n = $(letterPos[Alt] +1);     
    }
    printf("%s\t%s\t%s\t%s\n", depth_p, mis_p, depth_n, mis_n)
}'