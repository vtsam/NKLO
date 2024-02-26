#!/bin/bash

echo SPIN:
read N
echo NR OF TOTAL DERIVATIVES:
read k
echo NR OF GLUONS:
read K
echo INCLUDE PERMUTATIONS:
read perms

if [ $k -gt $[N-1] ]
then
	echo ERROR: Number of derivatives must be smaller than N-1

else

form -l -q -d N="$N" -d k="$k" -d K="$K" -d PERM="$perms" NKLO

sed -i 's/sum/sum_/g' NKLO.log 
sed -i 's/RuleNKLO/L RuleNKLO/g' NKLO.log

mv NKLO.log tmpNKLO.h

form -l -q -d N="$N" -d k="$k" -d K="$K" sums 

mv sums.log outputN"$N"k"$k"K"$K".h

fi
