#import "utils.typ": *



= Appendix A: Derivations <appendix-a>

The following derivation closely follows the presentation in #cite(label("cf")) and the reader is referred there for more details.

== Constructing $k_+(dot.c, dot.c)$ and $cal(H)_+$
Given a reproducing kernel Hilbert space (RKHS) $cal(H)$ with kernel $k(dot.c, dot.c)$, we construct a new RKHS $cal(H)_+$ with kernel $k_+(dot.c, dot.c)$ by applying the Stein operator to both arguments of the base kernel.
#enum(
  numbering: "Step 1.",
  [*$k(dot.c, dot.c) -> k_0(dot.c, dot.c).$* 
  Suppose that $phi in cal(H)$, then $S_p [phi] in cal(H)_0$, the RKHS with kernel
   $ k_0 (z, z') &= cal(S)_p^z cal(S)_p^(z') [k(z, z')], $
  where $cal(S)_p^z$ and $cal(S)_p^(z')$ denote the Stein operator acting on the first and second arguments of $k(dot.c, dot.c)$ respectively. This can be expressed explicitly as
  $ k_0 (z, z') = nabla_z dot.c nabla_(z') k(z, z') + u(x) dot.c nabla_(z') k(z, z') + u(z') dot.c nabla_z k(z, z') + u(z) dot.c u(z') k(z, z'). $
  It can be shown that $k_0(dot.c, dot.c)$ is a valid reproducing positive-definite kernel[cite] and  satisfies the _zero-mean property_
  $ integral_Omega k_0 (z, z') p(z') dif z' = 0 $
  for almost all $z in Omega$.
  Hence, any function $psi in cal(H)_0$ has zero expectations under $p$, i.e., $EE_p [psi(Z)] = 0$.],
  [*$cal(H)_0 -> cal(H)_+. $* The kernel $k_0(dot.c, dot.c)$ allows us to approximate target zero mean functions $h(z)$ with $psi in cal(H)_0$. To extend this to general $h(z)$, the mean component of $h(z)$ must be represented. 
  To this end, $cal(H)_0$ is modified with the RKHS of constant functions
  $ cal(C) = {c : c in RR}, $
  with the reproducing kernel $k_C (z, z') = 1$ for all $z, z' in Omega$.
  The new RKHS is defined as
  $ cal(H)_+ = cal(C) + cal(H)_0 = {c + psi : c in RR, psi in cal(H)_0} $
  and the corresponding kernel is simply the sum of the component kernels:
  $ k_+ (z, z') = k_C (z, z') + k_0 (z, z') = 1 + k_0 (z, z'). $
  ]
)
== Boundedness of the induced kernel $k_0$
The base kernel is the RBF kernel modified with an polynomial decay term.
$ k(z, z') = exp(-||z - z'||^2\/2 alpha_2^2) / (1 + alpha_1 (||z||^2 +||z'||^2)), quad  alpha_1, alpha_2 > 0. $
The density function $p_Z$ considered is the standard normal distribution on $RR^d$:
$ p(z) = (2 pi)^(-d/2) exp(- ||z||^2\/ 2). $

#enum(
  numbering: "Step 1.",
  [Compute the score function for $p_Z$:
  #align(center)[
  $ u(z) = nabla_Z dot.c log p(z) = -z => cal(S)_p^z [phi(z)] = nabla_z phi(z) - z dot.c phi (z). $]
  ],
  [Apply Stein operator to left argument $z$:
  $  k_0(z, z') = cal(S)_p^(z') cal(S)_p^(z) [k(z, z')] $
  $ phi(z) = k(z, dot.c) ==> nabla_z phi(z) &= -exp(-||z - z'||^2\/2 alpha_2^2) / (1 + alpha_1 (||z||^2 +||z'||^2)) ((z - z') / (alpha_2^2) - (2 alpha_1 z) / (1 + alpha_1( ||z||^2 + ||z'||^2))) \ &= -k(z, z') ((z - z') / (alpha_2^2) - (2 alpha_1 z) / (1 + alpha_1( ||z||^2 + ||z'||^2))). $
  Hence, we have
  $ cal(S)_p [k (z, z')] &= -k(z, z') ((z - z') / (alpha_2^2) - (2 alpha_1 z) / (1 + alpha_1( ||z||^2 + ||z'||^2))) - z dot.c k(z,z') \
  &= -k(z, z') ((z - z') / (alpha_2^2) - (2 alpha_1 z) / (1 + alpha_1( ||z||^2 + ||z'||^2)) + z). $],
  [Apply on the right argument $z'$:
  $ cal(S)_p^(z') [cal(S)_p^(z) [k(z, z')]] = nabla_(z') [cal(S)_p^(z) [k(z, z')]] - z' dot.c [cal(S)_p^(z) [k(z, z')]] $
  $ dots.c = k(z, z') ( (1 + ||z - z'||^2)/(alpha_2^2) + (||z - z'||^2)/(alpha_2^4) + (2 alpha_1 (||z - z'||^2 + 2 z^T z'))/(1 + alpha_1 (||z||^2 + ||z'||^2)) + z ^T z') $
  ]
)
It is easy to see that $k_0(z, z')$ is bounded for all $z, z' in RR^d$ since the RBF kernel decays exponentially as $||z||, ||z'|| -> oo$ while the rational function terms all have denomiators that dominate (or are of the same order as) the numerators. Hence, there exists some constant $C > 0$ such that $sup_(z, z' in RR^d) k_0(z, z') <= C$.

#pagebreak()
= Appendix B: Code
All of the code used to produce the results in this report is open source and publicly available on my  #link("https://https://github.com/chinzhening/capstone-project")[#text(fill: blue.darken(10%))[github]]. The repository is structured as follows:
#list(

  [*`src/`*: contains all source code including Jupyter notebooks and utility scripts
  #list(
    [A stable implementation of the control functionals method in `numpy` is provided in `src/control_functional_implementation.py`.],
    [The examples from #cite(label("textbook")) are reproduced in the notebooks `ch_4.ipynb`, `ch_6.ipynb`, `ch_7.ipynb` and `ch_9.ipynb`. And selected numerical exercises are implemented in the notebooks titled `ex_[*].ipynb`, where `*` is the exercise number.],
    [For report specific experiments, see
    - `plain_monte_carlo.ipynb`
    - `antithetic_sampling.ipynb`
    - `control_variates.ipynb`
    - `control_functionals_rate_of_convergence.ipynb`
    - `control_functionals.ipynb`
    - `importance_sampling.ipynb`,
    -`portfolio_credit_risk.ipynb`]
  )
  ],
  [*`report/`*: contains the report source files written in Typst including bibliography and figures, as well as the compiled PDF version of the report and presentation slides.],
)
