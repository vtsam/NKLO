Auto S i,j,x,m,l,del,DP,Dp;
CF del,binom,pow,C;
Auto F f;
F sum;

#if 1==0

Implement the NKLO Feynman rule for the leading-twist quark operators.

There are 3 input parameters: 

- N: Spin of the operator
- k: Number of total derivatives
- K: Number of gluons

File by S. Van Thurenhout

#endif

Off Statistics;

#procedure iterSum

id sum(?a)*f(?b) = f(sum(?a,?b));

SplitArg f;

repeat id f(x1?,x2?,?a) = f(x1)+f(x2,?a);

#endprocedure

#if `K'=0

L RuleNKLO = (del2*DP2)^(`N'-`k'-1)*(del1*DP1+del2*DP2)^`k';

#else

#if `PERM'=0

L RuleNKLO = pow(-delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i`K',0,i{`K'-1}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m`K',0,i{`K'-1}-i`K'+m{`K'-1}-1)>*f(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-1}-i`K'+m{`K'-1}-1,m`K')>*<pow(Dp1,(i0-i1+m0-m1-1))>*...*<pow(Dp`K',(i{`K'-1}-i`K'+m{`K'-1}-m`K'-1))>*pow(del2*DP2,(`N'-`k'-i1+m`K'-2))*binom(`k',l)*pow(del1*DP1,`k'-l));

Transform,f,replace(1,last)=(binom,binom_);

Argument sum,f;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

Argument pow,binom_;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

EndArgument;

EndArgument;

#do i=1,2*`K'+1
#call iterSum
#enddo

id f(x?) = x;

#else

L RuleNKLO = perm_(fperm,<Dp1>,...,<Dp`K'>);

id fperm(?a) = pow(-delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i`K',0,i{`K'-1}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m`K',0,i{`K'-1}-i`K'+m{`K'-1}-1)>*f(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-1}-i`K'+m{`K'-1}-1,m`K')>*<pow(Dp1,(i0-i1+m0-m1-1))>*...*<pow(Dp`K',(i{`K'-1}-i`K'+m{`K'-1}-m`K'-1))>*pow(del2*DP2,(`N'-`k'-i1+m`K'-2))*binom(`k',l)*pow(del1*DP1,`k'-l))*fperm(?a);

Transform,f,replace(1,last)=(binom,binom_);

Argument sum,f;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

Argument pow,binom_;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

EndArgument;

EndArgument;

id fperm(x?) = fcol(Dp1);
id fperm(<x1?>,...,<x`K'?>) = replace_(<Dp1,x2>,...,<Dp{`K'-1},x`K'>,Dp`K',x1)*fcol(<Dp1>,...,<Dp`K'>);

#do i=1,2*`K'+1
#call iterSum
#enddo

id f(x?) = x;

#endif
#endif

P +f;
.end
