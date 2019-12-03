#!/bin/sh

mawk '
BEGIN {
    lines = 0;
    # naming for the strandedness info
    DPF[1]="R1+";
    DPF[2]="R1-";
    DPF[3]="R2+";
    DPF[4]="R2-";

    # ALL FIELD SWITCH #
    ALL = 0;
    ################ VCF2VARSCAN TRANSLATOR ###################
    # only the mentioned fields are used in the output 
    VCF["SS"] = "somatic_status";
    VCF["GPV"] = "variantP";
    VCF["SPV"] = "somaticP";
    VCF["N-RD"] = "NR1";
    VCF["N-AD"] = "NR2";
    VCF["T-RD"] = "TR1";
    VCF["T-AD"] = "TR2";
    VCF["N-R1+"] = "NR1+";
    VCF["N-R2+"] = "NR2+";
    VCF["T-R1+"] = "TR1+";
    VCF["T-R2+"] = "TR2+";

    VCF["somatic_status"] = "SS";
    VCF["variantP"] = "GPV";
    VCF["somaticP"] = "SPV";
    VCF["NR1"] = "N-RD";
    VCF["NR2"] = "N-AD";
    VCF["TR1"] = "T-RD";
    VCF["TR2"] = "T-AD";
    VCF["NR1+"] = "N-R1+";
    VCF["NR2+"] = "N-R2+";
    VCF["TR1+"] = "T-R1+";
    VCF["TR2+"] = "T-R2+";

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


    ######## HEADER OUTPUT #############
    printf("Chr\tStart\tEnd\tRef\tAlt");
    # all output 
    if (ALL == 1) {
        for (i in Data) {
            printf("%s\t",i);
        }
    # SELECT output
    } else {
        for (i = 4; i++ < 15;) {
            printf("\t%s",HEADER[i]);
        }
    }
    printf("\n");
}



##### PROCESS LINES #################
/^[^#]/ {
    lines++;
    start=$2;
    R=$4
    A=$5
    RL=length($4);
    AL=length($5);

    # convert VCF tags to array Data
    Info=$8;
    In = split(Info,IF,";");
    for (i = 0; ++i <= In;) {
        split(IF[i],s,"=");
        # convert the SS field value 1,2,3 to Germline, Somatic or LOH
        if (s[1] == "SS") {
            if (s[2] == 1) {
                Data[s[1]] = "Germline";
            } else if (s[2] == 2) {
                Data[s[1]] = "Somatic";
            } else {
                Data[s[1]] = "LOH";
            }
        } else {
            Data[s[1]] = s[2];
        }
        
    }
    # IF is Info tag array
    # normal values
    Format=$9;
    n = split(Format,FFields,":");
    N=$10;   
    split(N,NFields,":");
    # tumor values
    T=$11;
    split(T,TFields,":");
    # get the data values into an array
    # N-DP4 value
    # T-DP4 value etc.
    for (i = 0; ++i <= n;) {
        if (FFields[i] == "DP4") {
            split(TFields[i],TDP,",");
            split(NFields[i],NDP,",");
            for (j=0; j++ < 5;) {
                Data["T-" DPF[j]] = TDP[j];
                Data["N-" DPF[j]] = NDP[j];
            }
        } else {
            Data["T-" FFields[i]]=TFields[i];
            Data["N-" FFields[i]]=NFields[i]; 
        }

    }
    # ######### DEBUG THE DATA FIELDS ########
    # for (i in Data) {
    #     printf("%s:%s;",i,Data[i])
    # }
    # print "\n";
    # ########################################

    ######### CHANGE INDEL COORDS ############
    # variant is an insertion   
    if (AL > 1) {
        if (RL == 1) {
            end=start;
            # remove the redundant base in Ref and Alt
            R="-";
            A=substr(A,2)
        }
    # variant is a deletion
    } else if (RL > 1) {
        if (AL == 1) {
            # make the pos of deleted bases to start and end
            start=start+1;
            end=start + RL -2;
            # remove the redundant base in Ref and Alt
            R=substr(R,2);
            A="-";
        }
    } else {
        end = start;
    }

    ######## OUTPUT #################
    printf("%s\t%s\t%s\t%s\t%s",$1,start,end,R,A);
    if (ALL == 1) {
        for (i in Data) {
            printf("\t%s",Data[i]);
        }
    } else {
        for (i = 4; i++ < 15;) {
        printf("\t%s",Data[VCF[HEADER[i]]]);
        }
    }
    printf("\n");

}'