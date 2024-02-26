FORM implementation of the NKLO Feynman rule. There are 3 main files:

* NKLO.frm: Sets up the Feynman rule, leaving the sums unevaluated.
* sums.frm: Evaluates the sums
* NKLO.sh: Shell script to generate the Feynman rule. When run, the user will be prompted to input the operator spin, the number of total derivatives and the number of gluons. This will generate the kinematic part of the Feynman rule corresponding
  to the color-ordering t^{c_1}...t_{c_K}. The FORM output does NOT include the other permutations of the gluon fields, and also the overall factor g_s (t^{c_1}...t^{c_K})(Delta_{mu_1}...Delta_{mu_K})(Delta.Gamma) is omitted. The rule is written to
  the output file outputN...k...K... .h. For example, for N=3, k=1 and K=1 the output file is called outputN3k1K1.h and reads

  RuleNKLO =
       - Dp(1)*delD - DP(1)*del1*delD - DP(2)*del2*delD;

The notation is as follows:

* delD sets the conventions used for the sign of the strong coupling in the covariant derivative
* del1 and del2 set the conventions used for the momenta of the quark momenta, denoted by DP(1) and DP(2) resp. In particular, for all momenta flowing into the operator vertex, one should set del1=del2=+1. For physical momentum flow, set
  del1=-1, del2=+1.
* Dp(i) denotes the momentum of gluon i.

More details on the notation and conventions can be found in Sec.2 of the main text.
