Auto S a,b,c,i,j,x,m,l,del,DP,Dp,mu;
S gs;
CF del,binom,pow,C,Delta,DeltaGamma,dd,g;
Auto CF P,p,eps;
Auto F f;
F sum,t;

#if 1==0

Implement the NKLO Feynman rule for the leading-twist polarized gluon operators.

There are 3 input parameters: 

- N: Spin of the operator
- k: Number of total derivatives
- K: Number of additional gluons (the total nr is K+2)

File by S. Van Thurenhout

#endif

Off Statistics;

#procedure iterSum(f)

id sum(?a)*`f'(?b) = `f'(sum(?a,?b));

SplitArg `f';

repeat id `f'(x1?,x2?,?a) = `f'(x1)+`f'(x2,?a);

#endprocedure

#if `K'=0

#if `PERM'=0

L RuleNKLOgp = (-i_)^`N'*dd(a,b)*Delta(P2)^(`N'-`k'-2)*(Delta(P1)+Delta(P2))^`k'*(Delta(P2)*eps(Delta,mu1,mu2,P1)+Delta(mu2)*eps(Delta,mu1,P1,P2));

#else

L RuleNKLOgp = (-i_)^`N'*dd(a,b)*(Delta(P1)+Delta(P2))^`k'*(Delta(P2)^(`N'-`k'-2)*(Delta(P2)*eps(Delta,mu1,mu2,P1)+Delta(mu2)*eps(Delta,mu1,P1,P2))-Delta(P1)^(`N'-`k'-2)*(Delta(P1)*eps(Delta,mu1,mu2,P2)+Delta(mu1)*eps(Delta,mu2,P1,P2)));

#endif

#else

#if `PERM'=0

* double Abelian

L AbelRule = (Delta(P2)*eps(mu{`K'+2},Delta,P1,mu{`K'+1})-Delta(mu{`K'+2})*eps(P2,Delta,P1,mu{`K'+1}))*pow(gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i`K',0,i{`K'-1}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m`K',0,i{`K'-1}-i`K'+m{`K'-1}-1)>*fAbel(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-1}-i`K'+m{`K'-1}-1,m`K')>*<pow(-i_*Delta(p1),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p`K'),(i{`K'-1}-i`K'+m{`K'-1}-m`K'-1))>*pow(-i_*Delta(P2),(`N'-`k'-i1+m`K'-3))*binom(`k',l)*pow(-i_*Delta(P1),`k'-l))*<Delta(mu1)>*...*<Delta(mu`K')>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>;

Transform,fAbel,replace(1,last)=(binom,binom_);

Argument sum,fAbel;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-2;

Argument pow,binom_;

id m0 = i`K'+l-i0+i1+1;
id i0 = `N'-`k'-2;

EndArgument;

EndArgument;

#do i=1,2*`K'+1
#call iterSum(fAbel)
#enddo

id fAbel(x?) = x;

.sort


* Abelian non-Abelian

L NAbelRule = 0 +

#if `K'=1

((-1)^`K'*i_/2*(eps(mu{`K'+1},Delta,mu{`K'+2},mu{`K'})*Delta(P1)-eps(P1,Delta,mu{`K'+2},mu{`K'})*Delta(mu{`K'+1}))*pow(gs*delD,`K')*sum(l,0,`k')*fNAbel(pow(-i_*Delta(P1),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P2)-i_*Delta(p`K'),`k'-l))*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)

+

((-1)^`K'*i_*eps(mu{`K'+1},Delta,P2,mu{`K'+2})*pow(gs*delD,`K')*sum(l,0,`k')*fNAbel(pow(-i_*Delta(p1)-i_*Delta(P1),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P2),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>);

#else

((-1)^`K'*i_/2*(eps(mu{`K'+1},Delta,mu{`K'+2},mu{`K'})*Delta(P1)-eps(P1,Delta,mu{`K'+2},mu{`K'})*Delta(mu{`K'+1}))*pow(gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i{`K'-1},0,i{`K'-2}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m{`K'-1},0,i{`K'-2}-i{`K'-1}+m{`K'-2}-1)>*fNAbel(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-2}-i{`K'-1}+m{`K'-2}-1,m{`K'-1})>*<pow(-i_*Delta(p{`K'-1}),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p1),(i{`K'-2}-i{`K'-1}+m{`K'-2}-m{`K'-1}-1))>*pow(-i_*Delta(P1),(`N'-`k'-i1+m{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2)-i_*Delta(p`K'),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'-1})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)

+

((-1)^`K'*i_*eps(mu{`K'+1},Delta,P2,mu{`K'+2})*pow(gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i{`K'-1},0,i{`K'-2}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m{`K'-1},0,i{`K'-2}-i{`K'-1}+m{`K'-2}-1)>*fNAbel(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-2}-i{`K'-1}+m{`K'-2}-1,m{`K'-1})>*<pow(-i_*Delta(p`K'),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p2),(i{`K'-2}-i{`K'-1}+m{`K'-2}-m{`K'-1}-1))>*pow(-i_*Delta(P1)-i_*Delta(p1),(`N'-`k'-i1+m{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>);

#endif

Transform,fNAbel,replace(1,last)=(binom,binom_);

Argument sum,fNAbel;

id m0 = i{`K'-1}+l-i0+i1+1;
id i0 = `N'-`k'-2;

Argument pow,binom_;

id m0 = i{`K'-1}+l-i0+i1+1;
id i0 = `N'-`k'-2;

EndArgument;

EndArgument;

#do i=1,2*`K'+1
#call iterSum(fNAbel)
#enddo

id fNAbel(x?) = x;

.sort

* double non-Abelian

L DNAbelRule = 0 +

#if `K'=1

fDNAbel(0);

#elseif `K'=2

-1/2*eps(mu{`K'+2},Delta,mu{`K'+1},mu1)*pow(gs*delD,`K')*sum(l,0,`k')*fDNAbel(pow(-i_*Delta(p`K')-i_*Delta(P2),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P1)-i_*Delta(p1),`k'-l))*Delta(mu2)*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>;

#else

(-(-1)^`K'/2*eps(mu{`K'+2},Delta,mu{`K'+1},mu1)*pow(gs*delD,`K')*<sum(i1,0,i0-1)>*...*<sum(i{`K'-2},0,i{`K'-3}-1)>*sum(l,0,`k')*<sum(m1,0,i0-i1+m0-1)>*...*<sum(m{`K'-2},0,i{`K'-3}-i{`K'-2}+m{`K'-3}-1)>*fDNAbel(<binom(i0-i1+m0-1,m1)>*...*<binom(i{`K'-3}-i{`K'-2}+m{`K'-3}-1,m{`K'-2})>*<pow(-i_*Delta(p2),(i0-i1+m0-m1-1))>*...*<pow(-i_*Delta(p{`K'-1}),(i{`K'-3}-i{`K'-2}+m{`K'-3}-m{`K'-2}-1))>*pow(-i_*Delta(p`K')-i_*Delta(P2),(`N'-`k'-i1+m{`K'-2}-3))*binom(`k',l)*pow(-i_*Delta(P1)-i_*Delta(p1),`k'-l))*<Delta(mu2)>*...*<Delta(mu`K')>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>);

#endif

Transform,fDNAbel,replace(1,last)=(binom,binom_);

Argument sum,fDNAbel;

id m0 = i{`K'-2}+l-i0+i1+1;
id i0 = `N'-`k'-2;

Argument pow,binom_;

id m0 = i{`K'-2}+l-i0+i1+1;
id i0 = `N'-`k'-2;

EndArgument;

EndArgument;

#do i=1,2*`K'+1
#call iterSum(fDNAbel)
#enddo

id fDNAbel(x?) = x;

.sort

L RuleNKLOgp = AbelRule + NAbelRule + DNAbelRule;

id eps(x1?,Delta,?a) = -eps(Delta,x1,?a);

repeat;

#do i=1,{`K'+2}
#do j=1,{`K'+2}

#if `i'>`j'
	id eps(Delta,?a,mu`i',?b,mu`j',?c) = -eps(Delta,?a,mu`j',?b,mu`i',?c);
	id eps(Delta,?a,P`i',?b,P`j',?c) = -eps(Delta,?a,P`j',?b,P`i',?c);
	id eps(Delta,?a,p`i',?b,p`j',?c) = -eps(Delta,?a,p`j',?b,p`i',?c);
#endif

id eps(Delta,?a,P`i',mu`j',?b) = -eps(Delta,?a,mu`j',P`i',?b);
id eps(Delta,?a,p`i',mu`j',?b) = -eps(Delta,?a,mu`j',p`i',?b);
id eps(Delta,?a,P`i',p`j',?b) = -eps(Delta,?a,p`j',P`i',?b);

#enddo
#enddo

endrepeat;

.sort
Drop AbelRule,NAbelRule,DNAbelRule;

#else

L RuleNKLOgp = perm_(fperm,fp(P1,mu{`K'+1},a),fp(P2,mu{`K'+2},b),<fp(p1,mu1,c1)>,...,<fp(p`K',mu`K',c`K')>);

id fperm(?a) = (Delta(P2)*eps(mu{`K'+2},Delta,P1,mu{`K'+1})-Delta(mu{`K'+2})*eps(P2,Delta,P1,mu{`K'+1}))*pow(gs*delD,`K')*<sum(iab1,0,iab0-1)>*...*<sum(iab`K',0,iab{`K'-1}-1)>*sum(l,0,`k')*<sum(mab1,0,iab0-iab1+mab0-1)>*...*<sum(mab`K',0,iab{`K'-1}-iab`K'+mab{`K'-1}-1)>*fAbel(<binom(iab0-iab1+mab0-1,mab1)>*...*<binom(iab{`K'-1}-iab`K'+mab{`K'-1}-1,mab`K')>*<pow(-i_*Delta(p1),(iab0-iab1+mab0-mab1-1))>*...*<pow(-i_*Delta(p`K'),(iab{`K'-1}-iab`K'+mab{`K'-1}-mab`K'-1))>*pow(-i_*Delta(P2),(`N'-`k'-iab1+mab`K'-3))*binom(`k',l)*pow(-i_*Delta(P1),`k'-l))*<Delta(mu1)>*...*<Delta(mu`K')>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>*fperm(?a)

+ 0

#if `K'=1

+ ((-1)^`K'*i_/2*(eps(mu{`K'+1},Delta,mu{`K'+2},mu{`K'})*Delta(P1)-eps(P1,Delta,mu{`K'+2},mu{`K'})*Delta(mu{`K'+1}))*pow(gs*delD,`K')*sum(l,0,`k')*fNAbel(pow(-i_*Delta(P1),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P2)-i_*Delta(p`K'),`k'-l))*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a)

+

((-1)^`K'*i_*eps(mu{`K'+1},Delta,P2,mu{`K'+2})*pow(gs*delD,`K')*sum(l,0,`k')*fNAbel(pow(-i_*Delta(p1)-i_*Delta(P1),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P2),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a)

+ fDNAbel(0)*fperm(?a);

#elseif `K'=2

+ ((-1)^`K'*i_/2*(eps(mu{`K'+1},Delta,mu{`K'+2},mu{`K'})*Delta(P1)-eps(P1,Delta,mu{`K'+2},mu{`K'})*Delta(mu{`K'+1}))*pow(gs*delD,`K')*<sum(inab1,0,inab0-1)>*...*<sum(inab{`K'-1},0,inab{`K'-2}-1)>*sum(l,0,`k')*<sum(mnab1,0,inab0-inab1+mnab0-1)>*...*<sum(mnab{`K'-1},0,inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1)>*fNAbel(<binom(inab0-inab1+mnab0-1,mnab1)>*...*<binom(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1,mnab{`K'-1})>*<pow(-i_*Delta(p{`K'-1}),(inab0-inab1+mnab0-mnab1-1))>*...*<pow(-i_*Delta(p1),(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-mnab{`K'-1}-1))>*pow(-i_*Delta(P1),(`N'-`k'-inab1+mnab{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2)-i_*Delta(p`K'),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'-1})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>*fperm(?a))

+

((-1)^`K'*i_*eps(mu{`K'+1},Delta,P2,mu{`K'+2})*pow(gs*delD,`K')*<sum(inab1,0,inab0-1)>*...*<sum(inab{`K'-1},0,inab{`K'-2}-1)>*sum(l,0,`k')*<sum(mnab1,0,inab0-inab1+mnab0-1)>*...*<sum(mnab{`K'-1},0,inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1)>*fNAbel(<binom(inab0-inab1+mnab0-1,mnab1)>*...*<binom(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1,mnab{`K'-1})>*<pow(-i_*Delta(p`K'),(inab0-inab1+mnab0-mnab1-1))>*...*<pow(-i_*Delta(p2),(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-mnab{`K'-1}-1))>*pow(-i_*Delta(P1)-i_*Delta(p1),(`N'-`k'-inab1+mnab{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a)

+

(-1/2*eps(mu{`K'+2},Delta,mu{`K'+1},mu1)*pow(gs*delD,`K')*sum(l,0,`k')*fDNAbel(pow(-i_*Delta(p`K')-i_*Delta(P2),(`N'-`k'-2+l))*binom(`k',l)*pow(-i_*Delta(P1)-i_*Delta(p1),`k'-l))*Delta(mu2)*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a);

#else

+ ((-1)^`K'*i_/2*(eps(mu{`K'+1},Delta,mu{`K'+2},mu{`K'})*Delta(P1)-eps(P1,Delta,mu{`K'+2},mu{`K'})*Delta(mu{`K'+1}))*pow(gs*delD,`K')*<sum(inab1,0,inab0-1)>*...*<sum(inab{`K'-1},0,inab{`K'-2}-1)>*sum(l,0,`k')*<sum(mnab1,0,inab0-inab1+mnab0-1)>*...*<sum(mnab{`K'-1},0,inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1)>*fNAbel(<binom(inab0-inab1+mnab0-1,mnab1)>*...*<binom(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1,mnab{`K'-1})>*<pow(-i_*Delta(p{`K'-1}),(inab0-inab1+mnab0-mnab1-1))>*...*<pow(-i_*Delta(p1),(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-mnab{`K'-1}-1))>*pow(-i_*Delta(P1),(`N'-`k'-inab1+mnab{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2)-i_*Delta(p`K'),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'-1})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>*fperm(?a))

+

((-1)^`K'*i_*eps(mu{`K'+1},Delta,P2,mu{`K'+2})*pow(gs*delD,`K')*<sum(inab1,0,inab0-1)>*...*<sum(inab{`K'-1},0,inab{`K'-2}-1)>*sum(l,0,`k')*<sum(mnab1,0,inab0-inab1+mnab0-1)>*...*<sum(mnab{`K'-1},0,inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1)>*fNAbel(<binom(inab0-inab1+mnab0-1,mnab1)>*...*<binom(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-1,mnab{`K'-1})>*<pow(-i_*Delta(p`K'),(inab0-inab1+mnab0-mnab1-1))>*...*<pow(-i_*Delta(p2),(inab{`K'-2}-inab{`K'-1}+mnab{`K'-2}-mnab{`K'-1}-1))>*pow(-i_*Delta(P1)-i_*Delta(p1),(`N'-`k'-inab1+mnab{`K'-1}-3))*binom(`k',l)*pow(-i_*Delta(P2),`k'-l))*<Delta(mu1)>*...*<Delta(mu{`K'})>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a)

+

(-(-1)^`K'/2*eps(mu{`K'+2},Delta,mu{`K'+1},mu1)*pow(gs*delD,`K')*<sum(idnab1,0,idnab0-1)>*...*<sum(idnab{`K'-2},0,idnab{`K'-3}-1)>*sum(l,0,`k')*<sum(mdnab1,0,idnab0-idnab1+mdnab0-1)>*...*<sum(mdnab{`K'-2},0,idnab{`K'-3}-idnab{`K'-2}+mdnab{`K'-3}-1)>*fDNAbel(<binom(idnab0-idnab1+mdnab0-1,mdnab1)>*...*<binom(idnab{`K'-3}-idnab{`K'-2}+mdnab{`K'-3}-1,mdnab{`K'-2})>*<pow(-i_*Delta(p2),(idnab0-idnab1+mdnab0-mdnab1-1))>*...*<pow(-i_*Delta(p{`K'-1}),(idnab{`K'-3}-idnab{`K'-2}+mdnab{`K'-3}-mdnab{`K'-2}-1))>*pow(-i_*Delta(p`K')-i_*Delta(P2),(`N'-`k'-idnab1+mdnab{`K'-2}-3))*binom(`k',l)*pow(-i_*Delta(P1)-i_*Delta(p1),`k'-l))*<Delta(mu2)>*...*<Delta(mu`K')>*<fcol(x0,c1,x1)>*...*<fcol(x{`K'-1},c`K',x`K')>)*fperm(?a);

#endif

Argument fcol;

id x0 = a;
id x`K' = b;

EndArgument;

Transform,fAbel,fNAbel,fDNAbel,replace(1,last)=(binom,binom_);

Argument fAbel,fNAbel,fDNAbel,sum,pow,binom_;

id mab0 = iab`K'+l-iab0+iab1+1;
id iab0 = `N'-`k'-2;

id mnab0 = inab{`K'-1}+l-inab0+inab1+1;
id inab0 = `N'-`k'-2;

id mdnab0 = idnab{`K'-2}+l-idnab0+idnab1+1;
id idnab0 = `N'-`k'-2;

Argument pow,binom_;

id mab0 = iab`K'+l-iab0+iab1+1;
id iab0 = `N'-`k'-2;

id mnab0 = inab{`K'-1}+l-inab0+inab1+1;
id inab0 = `N'-`k'-2;

id mdnab0 = idnab{`K'-2}+l-idnab0+idnab1+1;
id idnab0 = `N'-`k'-2;

EndArgument;

EndArgument;

id fperm(fp(PP1?,mupa?,b1?),fp(PP2?,mupb?,b2?),<fp(pp1?,mup1?,a1?)>,...,<fp(pp`K'?,mup`K'?,a`K'?)>) = replace_(P1,PP1,P2,PP2,<p1,pp1>,...,<p`K',pp`K'>,mu{`K'+1},mupa,mu{`K'+2},mupb,<mu1,mup1>,...,<mu`K',mup`K'>,a,b1,b,b2,<c1,a1>,...,<c`K',a`K'>);

Argument fcol;
	id a = c0;
	id b = c{`K'+1};
EndArgument;

repeat;

id eps(x1?,Delta,?a) = -eps(Delta,x1,?a);

#do i=0,`K'+2
#do j=0,`K'+2

#if `i'>`j'
	id fcol(c`i',c`j',c?) = -fcol(c`j',c`i',c);
	id fcol(c?,c`i',c`j') = -fcol(c,c`j',c`i');
	id eps(Delta,?a,mu`i',?b,mu`j',?c) = -eps(Delta,?a,mu`j',?b,mu`i',?c);
	id eps(Delta,?a,P`i',?b,P`j',?c) = -eps(Delta,?a,P`j',?b,P`i',?c);
	id eps(Delta,?a,p`i',?b,p`j',?c) = -eps(Delta,?a,p`j',?b,p`i',?c);
#endif

id eps(Delta,?a,P`i',mu`j',?b) = -eps(Delta,?a,mu`j',P`i',?b);
id eps(Delta,?a,p`i',mu`j',?b) = -eps(Delta,?a,mu`j',p`i',?b);
id eps(Delta,?a,P`i',p`j',?b) = -eps(Delta,?a,p`j',P`i',?b);

#enddo
#enddo

endrepeat;

#do i=1,{`K'+2}

#if `K'>2

id fcol(c?,cc1?,xx1?)*<fcol(ap0?,bp0?,cp0?)>*...*<fcol(ap{`K'-3}?,bp{`K'-3}?,cp{`K'-3}?)>*fcol(xx2?,c0,c`i') = (-1)^`K'*fcol(c0,c`i',x1)*<fcol(ap0,bp0,cp0)>*...*<fcol(ap{`K'-3},bp{`K'-3},cp{`K'-3})>*fcol(x{`K'-1},c,cc1);
id fcol(c?,c{`K'+1},xx1?)*<fcol(ap0?,bp0?,cp0?)>*...*<fcol(ap{`K'-3}?,bp{`K'-3}?,cp{`K'-3}?)>*fcol(xx2?,cc1?,c`i') = (-1)^`K'*fcol(cc1,c`i',x1)*<fcol(ap0,bp0,cp0)>*...*<fcol(ap{`K'-3},bp{`K'-3},cp{`K'-3})>*fcol(x{`K'-1},c,c{`K'+1});

#elseif `K'=2

id fcol(c?,cc1?,x1?)*fcol(x2?,c0,c`i') = (-1)^`K'*fcol(c0,c`i',x1)*fcol(x2,c,cc1);

#endif

#enddo

Argument fcol;
	id c0 = a;
	id c{`K'+1} = b;
EndArgument;

#do i=1,2*`K'+1
#call iterSum(fAbel)
#call iterSum(fNAbel)
#call iterSum(fDNAbel)
#enddo

id fAbel(x?) = x;
id fNAbel(x?) = x;
id fDNAbel(x?) = x;

#endif
#endif

P +f;
.end
