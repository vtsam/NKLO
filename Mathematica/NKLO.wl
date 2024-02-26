(* ::Package:: *)

Options[NKLO]={"kinematics"->"generic","covD"->"generic","PREF"->True};
NKLO[N_,k_,K_,opts:OptionsPattern[]]:=Block[{tab,tmp,rule0,rule,list,sum,i,j,m,l,DP,Dp,del1,del2,delD,gs,DiracGamma,perms,permsRl,t,Delta,mu,c},
OptionValue["kinematics"];
OptionValue["covD"];
OptionValue["PREF"];
tab=Join[Reverse[Table[list[m[j+1],0,i[j]-i[j+1]+m[j]-1],{j,0,K-1}]],Reverse[Table[list[i[j+1],0,i[j]-1],{j,0,K-1}]],{list[l,0,k]}]//Reverse;
rule0=tab/.{List->sum,list->List};
rule0=rule0 (-delD)^K Binomial[k,l](del1 DP[1])^(k-l)(del2 DP[2])^(N-k-i[1]+m[K]-2)Product[Binomial[i[j]-i[j+1]+m[j]-1,m[j+1]]Dp[j+1]^(i[j]-i[j+1]+m[j]-m[j+1]-1),{j,0,K-1}]/.{a_ sum[b__]:>sum[a,b]};
rule0=rule0//.{i[0]->N-k-1,m[0]->i[K]+l-i[0]+i[1]+1};
rule0=rule0/.{sum->Sum};
If[OptionValue["PREF"]==True,rule0=gs^K Dot[Delta,DiracGamma](NonCommutativeMultiply@@Table[t[Dp[i]],{i,1,K}])(Times@@Table[Delta[Dp[i]],{i,1,K}]) rule0;
perms=Permutations[Table[Dp[i],{i,1,K}]];
If[Length[perms[[1]]]>1,permsRl=Table[Table[Rule[perms[[j,i]],perms[[2,i]]],{i,1,Length[perms[[1]]]}],{j,1,Length[perms]}],permsRl={{Dp[1]->Dp[1]}}];
rule=0;
Do[rule=rule+(rule0/.permsRl[[jj]]),{jj,1,Length[permsRl]}];
rule=rule/.{Delta[Dp[a_]]:>Delta[mu[a]],t[Dp[a_]]:>t[c[a]]}/.NonCommutativeMultiply[a_]:>a,
rule=rule0;];
Which[OptionValue["covD"]=="generic",rule=rule,OptionValue["covD"]==1,rule=rule/.{delD->1},OptionValue["covD"]==-1,rule=rule/.{delD->-1}];
Which[OptionValue["kinematics"]=="incoming",rule/.{del1->1,del2->1},OptionValue["kinematics"]=="physical",rule/.{del1->-1,del2->1},OptionValue["kinematics"]=="generic",rule]
]
