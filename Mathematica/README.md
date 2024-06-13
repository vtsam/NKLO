Mathematica implementation of the NKLO Feynman rules. Details on the conventions can be found in Sec. 2 of the main text. The rules are implemented in NKLO.wl. The latter contains two functions, called NKLO and NKLOg, which take the arguments N, k and K (in that order). The first function computes the Feynman rules for the quark operators, while the second one computes the rules for the gluon operators. While N and k can be either numeric or symbolic, K has to be numeric (i.e. an integer larger than or equal to zero). Furthermore, the following options are accepted:

## NKLO

(1) "kinematics": choose the momentum routing, which corresponds to setting $\delta_1$ and $\delta_2$ to specific values. The default value is "generic", which keeps $\delta_1$ and $\delta_2$ symbolic. Other allowed values are "incoming", which sets $\delta_1=\delta_2=-1$ and "physical" for which $\delta_1=-1, \delta_2=+1$.

(2) "covD": specifies the convention used for the covariant derivative, which corresponds to setting a value for $\delta_D$. The default option is "generic", which keeps $\delta_D$ symbolic. Alternatively one can set "covD" $\rightarrow$ 1 or "covD" $\rightarrow$ -1.

(3) "perms": include all permutations of the gluons. When set to True, all $K!$ permutations are included. When set to False, the rule only keeps one particular ordering, corresponding to $t^{c_1}\dots t^{c_K}$. The default value is True.

## NKLOg

(1) "polarization": choose whether to generate the rules for the polarized or unpolarized gluon operators. The default value is False.

(2) "Abelian": choose whether to work in an Abelian or a non-Abelian theory. When set to True, only the Abelian part of the field strength is kept. The default value is False, meaning the generated rule is valid for a non-Abelian model like QCD. Note that, when setting this option to True, one needs to be careful with the interpretation of the prefactors, which are written specifically for QCD. This responsibility is put on the user.

(3) "covD": specifies the convention used for the covariant derivative, which corresponds to setting a value for $\delta_D$. The default option is "generic", which keeps $\delta_D$ symbolic. Alternatively one can set "covD" $\rightarrow$ 1 or "covD" $\rightarrow$ -1.

(4) "perms": include all permutations of the gluons. When set to True, all $(K+2)!$ permutations are included. Furthermore, the anti-symmetry properties of the structure constants are implemented to obtain a minimal set of color structures in the output. When set to False, the rule only keeps one particular ordering corresponding to $f^{a c_1 x_1}f^{x_1 c_1 x_2}\dots f^{x_{K-1}c_K b}$. The default value is True.
