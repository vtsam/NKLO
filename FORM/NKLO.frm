Auto S i,j,x,m,l,del,P,p;
S gs;
CF del,binom,pow,C,Delta,DeltaGamma,mu;
Auto F f;
F sum,t;

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

L RuleNKLO = DeltaGamma*(i_*del1*Delta(P1))^(`N'-`k'-1)*(i_*del1*Delta(P1)+i_*del2*Delta(P2))^`k';

#else

#if `PERM'=0

L RuleNKLO = pow(-i_*gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i`K',0,i{`K'-1}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m`K',0,i{`K'-1}-i`K'+m{`K'-1}-1)>*f(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-1}-i`K'+m{`K'-1}-1,m`K')>*<pow(-i_*Delta(p1),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p`K'),(i{`K'-1}-i`K'+m{`K'-1}-m`K'-1))>*pow(i_*del1*Delta(P1),(`N'-`k'-i1+m`K'-2))*binom(`k',l)*pow(i_*del2*Delta(P2),`k'-l))*<Delta(mu(p1))>*...*<Delta(mu(p`K'))>*DeltaGamma*<t(p1)>*...*<t(p`K')>;

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

L RuleNKLO = perm_(fperm,<p1>,...,<p`K'>);

id fperm(?a) = pow(-i_*gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i`K',0,i{`K'-1}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m`K',0,i{`K'-1}-i`K'+m{`K'-1}-1)>*f(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-1}-i`K'+m{`K'-1}-1,m`K')>*<pow(-i_*Delta(p1),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p`K'),(i{`K'-1}-i`K'+m{`K'-1}-m`K'-1))>*pow(i_*del1*Delta(P1),(`N'-`k'-i1+m`K'-2))*binom(`k',l)*pow(i_*del2*Delta(P2),`k'-l))*fperm(?a)*<Delta(mu(p1))>*...*<Delta(mu(p`K'))>*DeltaGamma*<t(p1)>*...*<t(p`K')>;

Transform,f,replace(1,last)=(binom,binom_);

Argument sum,f;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

Argument pow,binom_;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-1;

EndArgument;

EndArgument;

id fperm(x?) = fcol(p1);
id fperm(<x1?>,...,<x`K'?>) = replace_(<p1,x2>,...,<p{`K'-1},x`K'>,p`K',x1)*fcol(<p1>,...,<p`K'>);

id fcol(?a) = 1;

#do i=1,2*`K'+1
#call iterSum
#enddo

id f(x?) = x;

#endif
#endif

P +f;
.end
