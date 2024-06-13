(* ::Package:: *)

NKLO::usage = "NKLO[N,k,K] generates the Feynman rule for the spin-N quark operator with k total derivatives and K gluons. Note that we have the condition N-1 \[GreaterEqual] K+k. The notation follows the notation introduced
in arXiv:2403.12623 [hep-ph]:

* The momenta of the external quarks have momenta P[1] and P[2] while the K additional gluons have momenta p[1],...,p[K], color indices c[1],...,c[K] and Lorentz indices mu[1],...,mu[K].. 

* Delta represents the lightlike vector \[CapitalDelta], while the generators of the color group are denoted as t[c[i]].

The following options are available:

* kinematics: choose the momentum routing, which corresponds to setting \!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\) and \!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\) in the accompanying paper to specific values. The default option is \"generic\", which keeps \!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\) and \!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\) symbolic. Other valid
options are \"incoming\", which sets \!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\)=\!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\)=-1 and \"physical\" for which \!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\)=-1, \!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\)=+1.

* covD: specify the convention used for the sign of the coupling in the covariant derivative defined as \!\(\*SubscriptBox[\(D\), \(\[Mu]\)]\) = \!\(\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\) - \!\(\*SubscriptBox[\(\[Delta]\), \(D\)]\)(I \!\(\*SubscriptBox[\(g\), \(\(s\)\(\\\ \)\)]\)\!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\)). The default option is \"generic\", which keeps \!\(\*SubscriptBox[\(\[Delta]\), \(D\)]\) symbolic. Other valid options are 
covD -> 1 and covD -> -1.

* perms: include all permutations of the gluons. The default value is True. When set to False, the rule only contains 1 particular ordering corresponding to \!\(\*SuperscriptBox[\(t\), SubscriptBox[\(c\), \(1\)]]\)... \!\(\*SuperscriptBox[\(t\), SubscriptBox[\(c\), \(K\)]]\)."

NKLOg::usage = "NKLOg[N,k,K] generates the Feynman rule for the spin-N gluon operator with k total derivatives and K+2 gluons. Note that we have the condition N-2 \[GreaterEqual] K+k. The notation follows the notation in the article:

* The momenta of two of the gluons are denoted by P[1] and P[2] with color indices a and b and Lorentz indices mu[K+1] and mu[K+2]. The K additional gluons have momenta p[1],...,p[K], color indices c[1],...,c[K] and 
Lorentz indices mu[1],...,mu[K].

* Delta represents the lightlike vector \[CapitalDelta], while the structure constants are denoted as f[c[i],c[j],c[k]].

The following options are available:

* polarization: when set to True, the Feynman rule for the polarized gluon operator is generated. The default value is False.

* Abelian: choose whether to work in an Abelian or a non-Abelian theory. When set to True, only the Abelian part of the field strength is kept. The default value is False.

* covD: specify the convention used for the sign of the coupling in the covariant derivative defined as \!\(\*SubscriptBox[\(D\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(F\), \(\[Rho]\[Sigma]\)], \(a\)]\) = \!\(\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(F\), \(\[Rho]\[Sigma]\)], \(a\)]\) + \!\(\*SubscriptBox[\(\[Delta]\), \(D\)]\)(\!\(\*SubscriptBox[\(g\), \(s\)]\)\!\(\*SuperscriptBox[\(f\), \(abc\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(A\), \(\[Mu]\)], \(a\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(F\), \(\[Rho]\[Sigma]\)], \(c\)]\)) . The default option is \"generic\", which keeps \!\(\*SubscriptBox[\(\[Delta]\), \(D\)]\) symbolic. Other valid options are 
covD -> 1 and covD -> -1.

* perms: include all permutations of the gluons. The default value is True. When set to False, the rule only contains 1 particular ordering corresponding to \!\(\*SuperscriptBox[\(f\), \(\*SubscriptBox[\(ac\), \(1\)] \*SubscriptBox[\(x\), \(1\)]\)]\)...\!\(\*SuperscriptBox[\(f\), \(\*SubscriptBox[\(x\), \(K - 1\)] \*SubscriptBox[\(c\), \(K\)] b\)]\)." 

(*quarks*)

Options[NKLO]={"kinematics"->"generic","covD"->"generic","perms"->True};
NKLO[N_,k_,K_,opts:OptionsPattern[]]:=Block[{tab,tmp,rule0,rule,list,sum,i,j,m,l,(*DP,Dp,*)del1,del2,delD,gs,DiracGamma,perms,permsRl,t,Delta,mu,c,P,p,moms,cols,ins,momPerms,colPerms,inPerms,momRls,colRls,inRls},
OptionValue["kinematics"];
OptionValue["covD"];
OptionValue["perms"];
NKLO::SymK = "K should be an integer!";
NKLO::Kl = "The condition N-1 \[GreaterEqual] K+k is not fulfilled!";
If[Not[IntegerQ[K]==True],Message[NKLO::SymK] && Abort[]];
If[NumberQ[N] &&NumberQ[k] && NumberQ[K] && Not[N-1 >= K+k],Message[NKLO::Kl]; 0,
tab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-1}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-1}]],{list[l,0,k]}]//Reverse;
rule0=tab/.{List->sum,list->List};
rule0=rule0 (-I delD)^K Binomial[k,l](I del2 Delta[P[2]])^(k-l)(I del1 Delta[P[1]])^(N-k-i[1]+m[K]-2)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]](-I Delta[p[j+1]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-1}]/.{a_ sum[b__]:>sum[a,b]};
rule0=rule0//.{i[0]->N-k-1,m[0]->i[K]+l-i[0]+i[1]+1};
rule0=rule0/.{sum->Sum};
If[OptionValue["perms"]==True,

rule0=gs^K Dot[Delta,DiracGamma](NonCommutativeMultiply@@Table[t[c[i]],{i,1,K}])(Times@@Table[Delta[mu[i]],{i,1,K}]) rule0 /.NonCommutativeMultiply[]:>1/.{1**a_:>a};

moms=Table[p[ii],{ii,1,K}];
cols=Table[c[ii],{ii,1,K}];

momPerms=Permutations[moms];
colPerms=Permutations[cols];

momRls=Table[Table[moms[[ii]]->momPerms[[jj,ii]],{ii,1,Length[moms]}],{jj,1,Length[momPerms]}];
colRls=Table[Table[cols[[ii]]->colPerms[[jj,ii]],{ii,1,Length[cols]}],{jj,1,Length[colPerms]}];

rule=Plus@@Table[rule0/.momRls[[ii]]/.colRls[[ii]],{ii,1,Length[momRls]}];

,

rule=gs^K Dot[Delta,DiracGamma](NonCommutativeMultiply@@Table[t[c[i]],{i,1,K}])(Times@@Table[Delta[mu[i]],{i,1,K}])rule0/.NonCommutativeMultiply[a_]:>a/.NonCommutativeMultiply[]:>1/.\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"NonCommutativeMultiply", "[", 
RowBox[{"1", ",", "a__"}], "]"}], "\\[Rule]", 
RowBox[{"NonCommutativeMultiply", "[", "a", "]"}]}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\);];
Which[OptionValue["covD"]=="generic",rule=rule,OptionValue["covD"]==1,rule=rule/.{delD->1},OptionValue["covD"]==-1,rule=rule/.{delD->-1}];
Which[OptionValue["kinematics"]=="incoming",rule/.{del1->-1,del2->-1},OptionValue["kinematics"]=="physical",rule/.{del1->-1,del2->1},OptionValue["kinematics"]=="generic",rule]
]]

(*gluons*)

Options[NKLOg]={"Abelian"->False,"covD"->"generic","perms"->True,"polarization"->False};
NKLOg[N_,k_,K_,opts:OptionsPattern[]]:=Block[{Abtab,NAbtab,DNAbtab,tmp,Abrule0,NAbrule0,DNAbrule0,rule,list,sum,i,j,m,l,P,p,g,DP,Dp,del1,del2,delD,gs,perms,permsRl,t,Delta,mu,c,dd,a,b,f,x,moms,cols,ins,momPerms,colPerms,
inPerms,momRls,colRls,inRls,eps},
OptionValue["Abelian"];
OptionValue["covD"];
OptionValue["perms"];
OptionValue["polarization"];
NKLOg::SymK = "K should be an integer!";
NKLOg::Kl = "The condition N-2 \[GreaterEqual] K+k is not fulfilled!";
If[Not[IntegerQ[K]==True],Message[NKLOg::SymK] && Abort[]];
If[NumberQ[N] &&NumberQ[k] && NumberQ[K] && Not[N-2 >= K+k],Message[NKLOg::Kl]; 0,

If[OptionValue["polarization"]==False,

(*Unpolarized Abelian K-th order vertex*)

Abtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-1}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-1}]],{list[l,0,k]}]//Reverse;
Abrule0=Abtab/.{List->sum,list->List};
Abrule0=Abrule0 ( delD)^K Binomial[k,l](-I Delta[P[1]])^(k-l)(-I Delta[P[2]])^(N-k-i[1]+m[K]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[j+1]](-I Delta[p[j+1]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-1}]/.{a_ sum[b__]:>sum[a,b]};
Abrule0=Abrule0//.{i[0]->N-k-2,m[0]->i[K]+l-i[0]+i[1]+1};
Abrule0=Abrule0/.{sum->Sum};
Abrule0=-((P[1][mu[K+2]]Delta[mu[K+1]](Delta[P[2]])-P[1].P[2]Delta[mu[K+1]]Delta[mu[K+2]]-(Delta[P[1]])(Delta[P[2]])g[mu[K+1],mu[K+2]]+P[2][mu[K+1]]Delta[mu[K+2]](Delta[P[1]])))Abrule0;

(*Unpolarized non-Abelian (K-1)-th order vertex*)

If[K==0||OptionValue["Abelian"]==True,NAbrule0=0,
NAbtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-2}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-2}]],{list[l,0,k]}]//Reverse;
NAbrule0=NAbtab/.{List->sum,list->List};
NAbrule0=2(NAbrule0 ( delD)^K Delta[mu[1]](-I)(g[mu[K+1],mu[K+2]]Delta[P[2]]-P[2][mu[K+1]]Delta[mu[K+2]])Binomial[k,l](-I(Delta[P[1]]+Delta[p[1]]))^(k-l)(-I Delta[P[2]])^(N-k-i[1]+m[K-1]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[j+2]](-I Delta[p[j+2]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-2}]/.{a_ sum[b__]:>sum[a,b]});
NAbrule0=NAbrule0//.{i[0]->N-k-2,m[0]->i[K-1]+l-i[0]+i[1]+1};
NAbrule0=NAbrule0/.{sum->Sum};];

(*Unpolarized double non-Abelian (K-2)-th order vertex*)

If[K==0||K==1||OptionValue["Abelian"]==True,DNAbrule0=0,
DNAbtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-3}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-3}]],{list[l,0,k]}]//Reverse;
DNAbrule0=DNAbtab/.{List->sum,list->List};
DNAbrule0=DNAbrule0 (delD)^K g[mu[K+1],mu[K+2]]Delta[mu[1]]Delta[mu[K]]Binomial[k,l](-I(Delta[P[1]]+Delta[p[1]]))^(k-l)(-I(Delta[p[K]]+Delta[P[2]]))^(N-k-i[1]+m[K-2]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[j+2]](-I Delta[p[j+2]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-3}]/.{a_ sum[b__]:>sum[a,b]};
DNAbrule0=DNAbrule0//.{i[0]->N-k-2,m[0]->i[K-2]+l-i[0]+i[1]+1};
DNAbrule0=DNAbrule0/.{sum->Sum};];

,

(*Polarized Abelian K-th order vertex*)

Abtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-1}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-1}]],{list[l,0,k]}]//Reverse;
Abrule0=Abtab/.{List->sum,list->List};
Abrule0=Abrule0 (delD)^K Binomial[k,l](-I Delta[P[1]])^(k-l)(-I Delta[P[2]])^(N-k-i[1]+m[K]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[j+1]](-I Delta[p[j+1]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-1}]/.{a_ sum[b__]:>sum[a,b]};
Abrule0=Abrule0//.{i[0]->N-k-2,m[0]->i[K]+l-i[0]+i[1]+1};
Abrule0=Abrule0/.{sum->Sum};
Abrule0=(Delta[P[2]]eps[mu[K+2],Delta,P[1],mu[K+1]]-Delta[mu[K+2]]eps[P[2],Delta,P[1],mu[K+1]])Abrule0;

(*Polarized non-Abelian (K-1)-th order vertex*)

If[K==0||OptionValue["Abelian"]==True,NAbrule0=0,
NAbtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-2}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-2}]],{list[l,0,k]}]//Reverse;
NAbrule0=NAbtab/.{List->sum,list->List};
NAbrule0=(NAbrule0 ( delD)^K (-I)/2( eps[mu[K+1],Delta,mu[K+2],mu[K]]Delta[P[1]]- eps[P[1],Delta,mu[K+2],mu[K]]Delta[mu[K+1]])Binomial[k,l](-I(Delta[P[2]]+Delta[p[K]]))^(k-l)(-I Delta[P[1]])^(N-k-i[1]+m[K-1]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[K-j-1]](-I Delta[p[K-j-1]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-2}]/.{a_ sum[b__]:>sum[a,b]})
+(NAbrule0 ( delD)^K (-I)/2(2  eps[mu[K+1],Delta,P[2],mu[K+2]]Delta[mu[1]])Binomial[k,l](-I Delta[P[2]])^(k-l)(-I(Delta[P[1]]+Delta[p[1]]))^(N-k-i[1]+m[K-1]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[K-j]](-I Delta[p[K-j]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-2}]/.{a_ sum[b__]:>sum[a,b]});
NAbrule0=NAbrule0//.{i[0]->N-k-2,m[0]->i[K-1]+l-i[0]+i[1]+1};
NAbrule0=(-1)^(K+1)NAbrule0/.{sum->Sum}];

(*Polarized double non-Abelian (K-2)-th order vertex*)

If[K==0||K==1||OptionValue["Abelian"]==True,DNAbrule0=0,
DNAbtab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-3}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-3}]],{list[l,0,k]}]//Reverse;
DNAbrule0=DNAbtab/.{List->sum,list->List};
DNAbrule0=DNAbrule0 (delD)^K *1/2eps[mu[K+2],Delta,mu[K+1],mu[1]]Delta[mu[K]]Binomial[k,l](-I(Delta[P[1]]+Delta[p[1]]))^(k-l)(-I(Delta[p[K]]+Delta[P[2]]))^(N-k-i[1]+m[K-2]-3)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Delta[mu[j+2]](-I Delta[p[j+2]])^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-3}]/.{a_ sum[b__]:>sum[a,b]};
DNAbrule0=DNAbrule0//.{i[0]->N-k-2,m[0]->i[K-2]+l-i[0]+i[1]+1};
DNAbrule0=(-1)^(K+1)DNAbrule0/.{sum->Sum};];];

rule=Abrule0+NAbrule0+DNAbrule0;

If[K==0,rule=dd[a,b]rule,rule=gs^K(NonCommutativeMultiply@@Table[f[x[jj],c[jj+1],x[jj+1]],{jj,0,K-1}])rule/.{NonCommutativeMultiply[a_]->a}/.{x[0]->c[0],x[K]->c[K+1]}];

If[OptionValue["perms"]==False,rule=rule/.{eps[a_,b_,c_,d_]:> Signature[List[a,b,c,d]]eps[Sort[{a,b,c,d}]]},

moms=Join[{P[1],P[2]},Table[p[ii],{ii,1,K}]];
cols=Join[{c[0],c[K+1]},Table[c[ii],{ii,1,K}]];
ins=Join[{mu[K+1],mu[K+2]},Table[mu[ii],{ii,1,K}]];

momPerms=Permutations[moms];
colPerms=Permutations[cols];
inPerms=Permutations[ins];

momRls=Table[Table[moms[[ii]]->momPerms[[jj,ii]],{ii,1,Length[moms]}],{jj,1,Length[momPerms]}];
colRls=Table[Table[cols[[ii]]->colPerms[[jj,ii]],{ii,1,Length[cols]}],{jj,1,Length[colPerms]}];
inRls=Table[Table[ins[[ii]]->inPerms[[jj,ii]],{ii,1,Length[ins]}],{jj,1,Length[inPerms]}];

rule=Plus@@Table[rule/.momRls[[ii]]/.inRls[[ii]]/.colRls[[ii]],{ii,1,Length[momRls]}];
rule=rule/.{dd[c[c1_],c[c2_]]/;c1>c2:>dd[c[c2],c[c1]]}/.{f[c[c1_],c[c2_],x[c3_]]/;c1>c2:>-f[c[c2],c[c1],x[c3]]}/.{f[x[c3_],c[c1_],c[c2_]]/;c1>c2:>-f[x[c3],c[c2],c[c1]]}//.{NonCommutativeMultiply[a___,Times[-1,b_],d___]:>-NonCommutativeMultiply[a,b,d]};
rule=rule/.{NonCommutativeMultiply[f[c[c1_],c[c2_],x[c3_]],b___,f[x[c4_],c[c5_],c[c6_]]]/;c1>c5:>(-1)^K NonCommutativeMultiply[f[c[c5],c[c6],x[c3]],b,f[x[c4],c[c1],c[c2]]]};
rule=rule//.{f[c[c1_],c[c2_],c[c3_]]/;c2>c3:>-f[c[c1],c[c3],c[c2]],f[c[c1_],c[c2_],c[c3_]]/;c1>c2:>-f[c[c2],c[c1],c[c3]]};
rule=rule//.{g[mu[m1_],mu[m2_]]/;m1>m2:>g[mu[m2],mu[m1]]};
rule=rule/.{eps[a_,b_,c_,d_]:> Signature[List[a,b,c,d]]eps[Sort[{a,b,c,d}]]};
];

rule=rule/.{c[0]->a,c[K+1]->b};

Which[OptionValue["covD"]=="generic",rule=rule,OptionValue["covD"]==1,rule=rule/.{delD->1},OptionValue["covD"]==-1,rule=rule/.{delD->-1}]

]]
