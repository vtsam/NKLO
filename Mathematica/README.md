Mathematica implementation of the NKLO Feynman rule. The rule is implemented in NKLO.wl. The NKLO.wl file contains a single function, called NKLO, which takes the arguments N, k and K (in that order). 
While N and k can be either numeric or symbolic, K has to be numeric (i.e. an integer larger than or equal to zero). Furthermore, the following options are accepted:

(1) "kinematics": choose the momentum routing, which corresponds to setting $\delta_1$ and $\delta_2$ to specific values. The default option is "generic", which keeps $\delta_1$ and $\delta_2$ symbolic. The other options are "incoming", which sets $\delta_1=\delta_2=+1$ and "physical" for which $\delta_1=-1, \delta_2=+1$.

(2) "covD": specifies the convention used for the covariant derivative, which corresponds to setting a value for $\delta_D$. The default option is "generic", which keeps $\delta_D$ symbolic. Alternatively one can set covD $\rightarrow 1$ or covD $\rightarrow -1$.

(3) "PREF": include the appropriate factors of
    $g_s^{K}\left(\prod_{j}t^{c_j}\Delta_{\mu_j}\right)(\Delta\cdot\Gamma)$
    taking into account all possible permutations of the gluonic quantities. The default option is "PREF" $\rightarrow$ True. When set to False, the prefactor is omitted and only one color-ordering, corresponding to $t^{c_1}\dots t^{c_K}$, is generated.
Within the package, the momenta are denoted by $DP[1]$ and $DP[2]$ for the quark momenta $\Delta\cdot P_1$ and $\Delta\cdot P_1$ and $Dp[i]$ for the gluonic ones $\Delta\cdot p_i$.

More details on the conventions can be found in Sec.2 of the main text.
