FORM implementation of the NKLO Feynman rule. There are 3 main files:

* NKLO.frm: Sets up the Feynman rule, leaving the sums unevaluated.
* sums.frm: Evaluates the sums
* NKLO.sh: Shell script to generate the Feynman rule. When run, the user will be prompted to input the operator spin, the number of total derivatives and the number of gluons. Furthermore, the user is asked whether or not to include permutations. When the latter is given a value 0, only the kinematic part of the Feynman rule corresponding
  to the color-ordering $t^{c_1}...t^{c_K}$ is generated. The FORM output does NOT include the other permutations of the gluon fields, and also the overall factor $g_s (t^{c_1}...t^{c_K})(\Delta_{\mu_1}\dots\Delta_{\mu_K})(\Delta\cdot\Gamma)$ is omitted. If the value for the permutations is different from 0, all permutations of the gluonic fields are generated, together with the proper prefactor. The rule is written to
  the output file outputN...k...K... .h. For example, for N=4, k=1 and K=2 the output file is called outputN4k1K2.h and reads

  RuleNKLO =
       + delD^2 * (
          + Dp1
          + Dp2
          + DP2 * del2
          + DP1 * del1
          );

  if the permutations are turned off and

  RuleNKLO =
        t(1) * t(2) * Delta(mu(1)) * Delta(mu(2)) * DeltaGamma * gs^2 * delD^2 * (
          + Dp1
          + Dp2
          + DP2 * del2
          + DP1 * del1
          )
       + t(2) * t(1) * Delta(mu(1)) * Delta(mu(2)) * DeltaGamma * gs^2 * delD^2 * (
          + Dp1
          + Dp2
          + DP2 * del2
          + DP1 * del1
          );

otherwise. The notation is as follows:

* delD sets the conventions used for the sign of the strong coupling in the covariant derivative
* del1 and del2 set the conventions used for the momenta of the quark momenta, denoted by DP1 and DP2 resp. In particular, for all momenta flowing into the operator vertex, one should set del1=del2=+1. For physical momentum flow, set
  del1=-1, del2=+1.
* Associated to gluon i we have the momentum dpi and a color generator t(i).
* Delta(mu(i)) is the lightlike vector with Lorentz index mu(i).
* DeltaGamma represents the dot product of Delta with some Dirac gamma matrix.

More details on the notation and conventions can be found in Sec.2 of the main text.
