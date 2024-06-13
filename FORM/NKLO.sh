#!/bin/bash

# File by S. Van Thurenhout

rm -f NKLO.log
rm -f tmpNKLO.h

echo OPERATOR:
read op
echo SPIN:
read N
echo NR OF TOTAL DERIVATIVES:
read k
echo NR OF ADDITIONAL GLUONS:
read K
echo INCLUDE PERMUTATIONS:
read perms

if [ $op = q ]
then

if [ $[k+K] -gt $[N-1] ]
then
	echo ERROR: Sum of number of derivatives and number of gluons must be smaller than N-1

else

form -l -q -d N="$N" -d k="$k" -d K="$K" -d PERM="$perms" NKLO

sed 's/sum/sum_/g' NKLO.log | sed 's/RuleNKLO/L RuleNKLO/g' > tmpNKLO.h

rm -f NKLO.log

form -l -q -d N="$N" -d k="$k" -d K="$K" sums 

rm -f tmpNKLO.h

mv sums.log output_"$op"op_N"$N"_k"$k"_K"$K".h

fi

elif [ $op = g ]
then

if [ $[k+K] -gt $[N-2] ]
then
	echo ERROR: Sum of number of derivatives and number of extra gluons must be smaller than N-2

else

form -l -q -d N="$N" -d k="$k" -d K="$K" -d PERM="$perms" NKLOg

sed 's/sum/sum_/g' NKLOg.log | sed 's/RuleNKLOg/L RuleNKLOg/g' > tmpNKLO.h

rm -f NKLOg.log

form -l -q -d N="$N" -d k="$k" -d K="$K" sums 

rm -f tmpNKLO.h

mv sums.log output_"$op"op_N"$N"_k"$k"_K"$K".h

fi

elif [ $op = gp ]
then

if [ $[k+K] -gt $[N-2] ]
then
	echo ERROR: Sum of number of derivatives and number of extra gluons must be smaller than N-2

else

form -l -q -d N="$N" -d k="$k" -d K="$K" -d PERM="$perms" NKLOgp

sed 's/sum/sum_/g' NKLOgp.log | sed 's/RuleNKLOgp/L RuleNKLOgp/g' > tmpNKLO.h

rm -f NKLOg.log

form -l -q -d N="$N" -d k="$k" -d K="$K" sums 

rm -f tmpNKLO.h

mv sums.log output_"$op"op_N"$N"_k"$k"_K"$K".h

fi
else

echo ERROR: Unknown operator type

fi
