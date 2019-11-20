#!/bin/sh
sed -i -r '
1 s/Otherinfo/NR1\tNR2\tNVAF\tN_geno\tTR1\tTR2\tTVAF\tT_geno\tsomatic_status\tvariantP\tsomaticP\tTR1+\tTR1-\tTR2+\tTR2-\tNR1+\tNR1-\tNR2+\tNR2-/
' $1
