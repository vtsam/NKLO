Auto S i,j,x,m,l,del;
CF del,binom,pow,Delta,mu,DeltaGamma;
Auto CF DP,Dp;
Auto F f;
F sum,t;
S gs;

Off Statistics;

#include tmpNKLO.h

id pow(x1?,x2?) = x1^x2;
id binom(x1?,x2?) = binom_(x1,x2);

id fcol(<Dp1?>,...,<Dp`K'?>) = gs^`K'*<Delta(mu(Dp1))>*...*<Delta(mu(Dp`K'))>*DeltaGamma*<t(Dp1)>*...*<t(Dp`K')>;

#do i=1,`K'
	id Delta(mu(Dp`i')) = Delta(mu(`i'));
	id t(Dp`i') = t(`i');
#enddo

B t,gs,delD,DeltaGamma,Delta;

P +f +s;
.end
