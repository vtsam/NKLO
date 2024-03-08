#!/bin/bash

# File by S. Van Thurenhout

rm -f NKLO.log
rm -f tmpNKLO.h

echo SPIN:
read N
echo NR OF TOTAL DERIVATIVES:
read k
echo NR OF GLUONS:
read K
echo INCLUDE PERMUTATIONS:
read perms

if [ $[k+K] -gt $[N-1] ]
then
	echo ERROR: Sum of number of derivatives and number of gluons must be smaller than N-1

else

form -l -q -d N="$N" -d k="$k" -d K="$K" -d PERM="$perms" NKLO

sed 's/sum/sum_/g' NKLO.log | sed 's/RuleNKLO/L RuleNKLO/g' > tmpNKLO.h

rm -f NKLO.log

form -l -q -d N="$N" -d k="$k" -d K="$K" sums 

rm -f tmpNKLO.h

mv sums.log outputN"$N"k"$k"K"$K".h

cat sums.log

fi
