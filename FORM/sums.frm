Auto S i,j,x,m,l,del,a,b,c;
CF del,binom,pow,Delta,mu,DeltaGamma,dd,g,eps;
Auto CF DP,Dp,P,p;
Auto F f;
F sum,t;
S gs;

Off Statistics;

#if 1==0

Evaluate sums

File by S. Van Thurenhout

#endif

#include tmpNKLO.h

id pow(x1?,x2?) = x1^x2;
id binom(x1?,x2?) = binom_(x1,x2);

#do i=1,`K'
	id Delta(mu(p`i')) = Delta(mu(`i'));
	id t(p`i') = t(c`i');
#enddo

Argument fcol;
	id x0 = a;
	id x`K' = b;
EndArgument;

B t,gs,delD,DeltaGamma,Delta,i_,del1,del2,fcol,g,eps;

P +f;
.end
