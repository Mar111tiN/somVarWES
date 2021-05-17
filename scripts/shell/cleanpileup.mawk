#!/bin/sh

#v1.1.0
### cleans samtools mpileup output
# first any end-of-read markers like ^I or ^[ will be removed
# cleanpileup takes extra care not to change any of the quality fields 
# .. which might accidentally contain such traces
# next, the A+12T for inserts and the T-12A for deletions will be converted to I/i and D/d

# creates headers
####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # output depth
        -d|--output_depth)
        outputDepth=1;
        shift;
        ;;
        # reduce output samples
        # required if the input is mpileup for tumor and normal (only normal is required)
        -s|--use-samples)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            useSamples=$2;
            shift 2;
        else
            echo "<cleanpileup> Error: sample count argument is missing\n[-s|--use-samples (default=all)]" >&2
            exit 1;
        fi;;
        -*|--*=) # unsupported flags
        echo "<cleanpileup> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done



mawk '
NR == 1 {
    printDepth="'$outputDepth'";
    samples = (NF - 3) / 3;
    useSamples="'${useSamples:-0}'";
    if (useSamples > 0) {
        samples = (useSamples > samples) ? samples : useSamples;
    }
    printf("Chr\tStart\tRef");
    for (s=0;s++<samples;) {
        if (printDepth) {
            printf("\tDepth%s\tRead%s",s,s);
        } else {
            printf("\tRead%s",s);
        }
        
    }
    printf("\n");
}
{   
    # print generic fields
    printf($1);
    for (c=1;c++<3;) {
        printf("\t%s",$c);
    }

    # cycle through the reads
    for (i = 0; i++ < samples;) {
        col = (i*3) + 2;
        # remove position traces from all read fields
        gsub(/\^[^\t]|\$/,"",$col);
        # 
        while (match($col,/[+-][0-9]+/)) {
            pre = substr($col,1,RSTART-2);
            indel = substr($col,RSTART,1); # + or -
            base = substr($col,RSTART-1,1); # base before the deletion
            l = substr($col,RSTART+1,RLENGTH-1);
            indel_bases = substr($col,RSTART+RLENGTH,l);
            post = substr($col,RSTART+RLENGTH+l);
            if (indel == "-") {
                if (match(indel_bases,/[ACGT]/)) {
                    base = "D";
                } else {
                    base = "d";
                }
            } else {
                if (match(indel_bases,/[ACGT]/)) {
                    base = "I";
                } else {
                    base = "i";
                }            
            }
            $col = pre base post;
        }
        if (printDepth) {
            printf("\t%s\t%s",$(col-1),$col);
        } else {
            printf("\t%s",$col);
        }
    }
    printf("\n");
}   '