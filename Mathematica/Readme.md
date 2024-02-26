Mathematica implementation of the NKLO Feynman rule. The rule is implemented in NKLO.wl. The NKLO.wl file contains a single function, called NKLO, which takes the arguments N, k and K (in that order). Furthermore, it accepts the following options:
\begin{itemize}
    \item "\textit{kinematics}": choose the momentum routing, which corresponds to setting $\delta_1$ and $\delta_2$ above to specific values. The default option is "\textit{generic}", which keeps $\delta_1$ and $\delta_2$ symbolic. The other options are "\textit{incoming}", which sets $\delta_1=\delta_2=+1$ and "\textit{physical}" for which $\delta_1=-1, \delta_2=+1$.
    \item "\textit{covD}": specifies the convention used for the covariant derivative, which corresponds to setting a value for $\delta_D$ above. The default option is "\textit{generic}", which keeps $\delta_D$ symbolic. Alternatively one can set covD $\rightarrow 1$ or covD $\rightarrow -1$.
    \item "\textit{PREF}": include the appropriate factors of
    \begin{equation}
    \label{eq:PREF}
        g_s^{K}\left(\prod_{j=1}^{K}t^{c_j}\Delta_{\mu_j}\right)(\Delta\cdot\Gamma),
    \end{equation}
    taking into account all possible permutations of the gluonic quantities. The default option is\newline PREF $\rightarrow$ True. When set to False, the prefactor in Eq.(\ref{eq:PREF}) is omitted and only one color-ordering, corresponding to $t^{c_1}\dots t^{c_K}$, is generated.
\end{itemize}
Within the package, the momenta are denoted by $DP[1]$ and $DP[2]$ for the quark momenta $\Delta\cdot P_1$ and $\Delta\cdot P_1$ and $Dp[i]$ for the gluonic ones $\Delta\cdot p_i$.
