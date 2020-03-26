#!/bin/sh

tumor=$1;
pon=$2;
output=$3;

echo $tumor > ${output}.pre;
cat $pon >> ${output}.pre;

cat ${output}.pre | \
awk  '
BEGIN {
    FS="/";
}
# print tumor bam and store basename (before _A) in pattern
NR==1 {
    string=$NF;
    split(string,pat,"_");
    pattern="/" pat[1];
    print $0;
}
$0 !~ pattern {
    print $0;
}' > $output;
rm ${output}.pre;