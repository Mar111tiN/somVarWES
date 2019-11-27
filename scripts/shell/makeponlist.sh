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
NR==1 {
    string=$NF;
    split(string,pat,"_");
    print $0;
}
$0 !~ pat[1] {
    print $0;
}' > $output;
rm ${output}.pre;