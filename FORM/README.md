FORM implementation of the NKLO Feynman rules. The rules themselves are implemented in the following files:

* NKLO.frm: Sets up the Feynman rule for quark operators.
* NKLOg.frm: Sets up the Feynman rule for unpolarized gluon operators.
* NKLOgp.frm: Sets up the Feynman rule for polarized gluon operators.

In each of these, the sums that appear are left unevaluated, as the summation indices are symbolic in intermediate steps. They are then explicitly evaluated in the sums.frm file. 

The main file for the user is NKLO.sh, which is a shell script to generate the desired operator Feynman rule. When run, the user will be prompted to input the operator type, operator spin, the number of total derivatives and the number of additional gluons (i.e. the order of the strong coupling $g_s$). Valid values for the operator type are q, g and gp which respectively generate the Feynman rules for quark, unpolarized and polarized gluon operators. The values of the other options need to be numeric. Furthermore, the user is asked whether or not to include permutations. When the latter is given a value 0 the part of the Feynman rule corresponding to the color-ordering $t^{c_1}\dots t^{c_K}$ is generated for the quark operators and $f^{a c_1 x_1}f^{x_1 c_1 x_2}\dots f^{x_{K-1}c_K b}$ for the gluon operators. All other permutations are omitted. If the value for the permutations is different from 0, all permutations of the gluonic fields are generated. The rule is written to the output file output_...op_N..._k..._K... .h. For example, the quark rule for N=4, k=1 and K=2 the output file is called output_qop_N4_k1_K2.h and reads

RuleNKLO =

       + t(c1)*t(c2)*Delta(p1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c1)*t(c2)*Delta(p2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c1)*t(c2)*Delta(P1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del1*delD^2 * (  - 1 )

       + t(c1)*t(c2)*Delta(P2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del2*delD^2 * (  - 1 );

  if the permutations are turned off and

RuleNKLO =

       + t(c1)*t(c2)*Delta(p1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c1)*t(c2)*Delta(p2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c1)*t(c2)*Delta(P1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del1*delD^2 * (  - 1 )

       + t(c1)*t(c2)*Delta(P2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del2*delD^2 * (  - 1 )

       + t(c2)*t(c1)*Delta(p1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c2)*t(c1)*Delta(p2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      delD^2 * ( 1 )

       + t(c2)*t(c1)*Delta(P1)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del1*delD^2 * (  - 1 )

       + t(c2)*t(c1)*Delta(P2)*Delta(mu(1))*Delta(mu(2))*DeltaGamma*i_*gs^2*
      del2*delD^2 * (  - 1 );

otherwise. The notation is as follows:

* delD sets the conventions used for the sign of the strong coupling (denoted in the code by gs) in the covariant derivative.
* del1 and del2 set the conventions used for the routing of the quark momenta, denoted by P1 and P2 resp. In particular, for all momenta flowing into the operator vertex, one should set del1=del2=-1. For physical momentum flow, set
  del1=-1, del2=+1.
* For the quark operator rules, gluon j we have the Lorentz index mu(j), momentum pj and a color generator t(cj). 
* For the gluon operators, two of the gluon fields have momenta P1 and P2 with color indices a and b and Lorentz indices mu(K+1) and mu(K+2). The additional K gluons follow the same notation as for the quark operator rules. 
* The structure constants are denoted by fcol(a,b,c). dd(a,b) represents the Kronecker-delta in color space.
* Delta represents the arbitrary lightlike vector $\Delta$ (with $\Delta^2=0$).
* DeltaGamma represents ($\Delta\cdot\Gamma$).
* eps(a,b,c,d) represents the Levi-Civita symbol while g(mu(1),mu(2)) is the metric tensor.

More details on the notation and conventions can be found in Sec. 2 of the main text.
