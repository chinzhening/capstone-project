// Capstone Report
// Author: Chin Zhe Ning
// Date: 2025-10-21

#import "utils.typ": *

#import "@preview/algo:0.3.6": algo, i, d, comment, code
#import "@preview/lilaq:0.4.0" as lq
#import "@preview/zebraw:0.5.5": *

// setup
#set page(margin: (x: 2.5cm, y: 2.5cm), numbering: "1")
#set heading(numbering: "1.")
#set par(first-line-indent: 1em, spacing: 1em,)
#set enum(indent: 1.25em, spacing: 1.25em)

#set text(font: "CMU Serif")

#show heading: set text(
  size: 12pt,

  fill: aqua.darken(70%))
#show heading: it => {
  v(0.5em)
  it
}

#set enum(indent: 1.5em, spacing: 1.5em)


#cover_page(
  title: "Monte Carlo Methods and their Applications to Finance",
  student: "Chin Zhe Ning",
  supervisor: "Ren Weiqing",
  examiner: "",
  department: "Department of Mathematics",
  institution: "National University of Singapore",
  info: "Capstone Project for Semester 1, AY25/26",
)

#pagebreak()
#outline()

#pagebreak()

= Introduction

In modern finance, a derivative is a security whose value depends on the price of an underlying asset, a group of assets, or a benchmark.
Accurate pricing is essential, as even small mispricings can lead to significant financial consequences.
Derivatives are often highly sensitive to fluctuations in their underlying assets, which makes it difficult to value.
Prices are typically determined via risk-neutral valuation, in which the expected discounted payoff is computed under a risk-neutral measure.
For derivatives with simple payoffs, analytical solutions are well-established.
A canonical example is the Black-Scholes formula for a European call option, which provides a closed-form solution under the assumption that the underlying asset price follows a geometric Brownian motion.
However, for more complex derivatives, e.g., those that are path-dependent or involve multiple underlying assets, analytical solutions are either difficult or impossible to obtain.

In such cases, Monte Carlo methods provide a powerful framework for estimating expected payoffs.
These methods are particularly well-suited for path-dependent derivatives and derivatives involving multiple underlying assets, where conventional numerical techniques may be infeasible.
Despite their versatility, Monte Carlo methods exhibit certain limitations.
Their slow rate of convergence and the high computational cost required for precise estimates can be significant constraints.
Nevertheless, due to their robustness and broad applicability, Monte Carlo simulations have become the workhorse of computational finance.

This report provides an overview of Monte Carlo methods in finance. Working knowledge of derivative pricing is assumed
and the focus is on the implementation of Monte Carlo methods for pricing derivatives. We first introduce some basic concepts in option pricing, followed by the basic Monte Carlo method for estimating expected payoffs. The discussion then addresses the limitations of the basic method and also the standard variance reduction techniques.

#pagebreak()
= Option Pricing

Suppose that the payoff of a derivative is given by a function $h(X)$ of a random
variable $X$ related to the underlying asset. The price of the derivative is given
by the expected discounted payoff under a risk-neutral measure, i.e., $v = EE[e^(-r T) h(X)]$,
where $r$ is the risk-free interest rate and $T$ is the time to maturity.

The random variable $X$ can take various forms depending on the type of derivative. For example, for a European option, $X$ could be the stock price at maturity $S_T$; for an Asian option, $X$ stock prices at multiple time points; for a basket option, $X$ could be a vector of prices of multiple underlying assets. The underlying asset prices are often modeled using stochastic processes. A common model is the geometric Brownian motion
$  S_t = S_0 exp { (r - 1/2 sigma^2) t + sigma  W_t }, quad W_t ~ N(0, t) $
where $S_0$ is the initial stock price, $sigma$ is the volatility, and $W_t$ is a standard Brownian motion. Unless otherwise specified, we will assume that the underlying asset prices follow this model.

== Examples <option-examples>
#enum(
  numbering: "1.",
  [*European Call Option.* For a European call option with underlying stock price $S_t$ the payoff of the option at time $T$ is given by
  $ h(S_T) = (S_T - K)^+$,
  where $K$ is the strike price, and $z^+ = max(z, 0)$. The price of the option is
  $ v = EE[(S_T - K)^+]. $
  It is possible to solve for $v$ analytically from the Black-scholes formula[cite]. Given the parameters $S_0, K, r, sigma, T$, the price of the option is
  $ v = S_0 Phi(sigma sqrt(T) - theta) - K e^(-r T) Phi(-theta) $
  where 
  $ theta = 1 / (sigma sqrt(T)) log K / S_0 + (sigma / 2 - r / sigma)sqrt(T). $],
  [*Asian Call Option.* For a discretely monitored average price call option, the
  payoff is given by
  $ h(S_(t_1), S_(t_2), ..., S_(t_m)) = (1/m sum_(i= 1)^m S_(t_i) - K)^+, $
  where $0 < t_1 < dots.c < t_m = T$ are a fixed set of dates. Under additional assumptions,
  it is possible to derive a closed-form solution[cite], but in general, there is no analytical
  solution. In order to compute $v$, we can simulate discrete paths of the stock
  price ${S_(t_i)}_(i = 1)^m$, compute the discounted payoff
  and then average over $N$ simulated paths.],
)

To estimate the expected payoff $v$, we can use Monte Carlo methods, which rely on random sampling to approximate the expected value.

#pagebreak()
= Basic Monte Carlo
To estimate the quantity
$ v = EE[h(X)] $
where $h$ is a function of $X$, we generate $N$ independent and identically distributed (i.i.d.) samples ${X_i}_(i=1)^N$ from the distribution of $X$ with density $f_X$.
The plain Monte Carlo estimate of $v$ is given by
$ hat(v) = 1/ N sum_(i=1)^(N) h(X_i). $
This is an unbiased estimate $hat(v)$ of the true mean $v$, i.e., $EE[hat(v)] = v,$
and the variance of $hat(v)$ is given by
$ "Var"(hat(v)) = 1/ N^2 sum_(i=1)^N "Var"(h(X_i))  = "Var"(h(X)) / N, $
As the number of samples $N$ increases, the variance of $hat(v)$ decreases proportionally to $1\/N$. Consequently, the estimator 
$hat(v)$ converges to the true mean $v$ _in the long run_.

#algo(
  line-numbers: false,
  radius: 5pt,
  stroke: 0pt,
  fill: none,
  inset: (x: 3em, y: 1em),

  indent-size: 2em,
  header: [*Pseudocode:*],
  
)[
  for $i = 1, ..., N$:#i\
    generate $X_i$ from $f_X$ #d #i\
    set $H_i = h(X_i)$ #d\
  
  compute the estimate $ hat(v) = 1/N sum_(i = 1)^N H_i $\
  compute the standard error $ "S.E." = sqrt(1/(N(N-1)) (sum_(i = 1)^N X_i^2- N hat(v)^2)) $
]

== Examples
#enum(
  numbering: "1.",
  [*European Call Option.* Consider the problem of estimating the price of a European call option described in @option-examples with strike price $K = 60$ and maturity $T = 1$ year. Assuming the parameters
  $ S_0 = 50, r = 0.05, sigma = 0.2, $ we estimate the option price using the plain Monte Carlo method for various sample sizes $N$. The results are summarised in @plain-monte-carlo-results.

  #figure(
    caption: [European Call Option: Plain Monte Carlo Estimate vs Sample Size],
  )[
    #table(
      columns: 6,
      align: center,
      [Sample Size], [$N = 10^2$], [$N = 10^3$], [$N = 10^4$], [$N = 10^5$], [$N = 10^6$],
      [Estimate], [1.8824], [1.6428], [1.6586], [1.6119], [1.6181],
      [S.E.], [0.6072], [0.1388], [0.0447], [0.0135], [0.0043],
      [R.E.], [32.25%], [8.45%], [2.69%], [0.84%], [0.27%],
    )
  ] <plain-monte-carlo-results>

  As $N$ increases, the estimate $hat(v)$ converges to the true option price $v = 1.6237$, and the standard error decreases to zero, this demonstrates the consistency of the Monte Carlo estimate.],
  [*Asian Call Option.* Consider the problem of estimating the price of a discretely monitored Asian call option described in @option-examples. Assuming the parameters
  $ S_0 = 50, r = 0.05, sigma = 0.2, K = 60, T = 1 $
  and monitoring at $m = 100$ equally spaced time points, we estimate the option price using the plain Monte Carlo method for various sample sizes $N$. The results are summarised in @asian-monte-carlo-results.
  
  #figure(
    caption: [Asian Call Option: Plain Monte Carlo Estimate vs Sample Size],
  )[
    #table(
      columns: 6,
      align: center,
      [Sample Size], [$N = 10^2$], [$N = 10^3$], [$N = 10^4$], [$N = 10^5$], [$N = 10^6$],
      [Estimate], [0.2823], [0.2534], [0.2650], [0.2715], [0.2700],
      [S.E.], [0.1224], [0.0357], [0.0122], [0.0040], [0.0013],
      [R.E.], [43.36%], [14.09%], [4.60%], [1.47%], [0.47%],
    )
  ] <asian-monte-carlo-results>

  The standard error decreases as the sample size $N$ increases, as observed in the European call option example. 
  ]
)

== Limitations
The numerical results in the previous section for the European call and Asian options demonstrate that, although the plain Monte Carlo method provides unbiased estimates, its convergence rate is relatively slow.

According to the Central Limit Theorem (CLT), if $"Var"(h(X)) = sigma^2 < oo$, then as $N -> oo$,
$ sqrt(N) (hat(v) - v) -> N(0, sigma^2), $

implying that for large $N$, the estimation error $hat(v) - v$ is approximately normally distributed with mean 0 and standard deviation $sigma \/ sqrt(N)$. So to be precise, the standard error (S.E.) of the estimate decays at rate $cal(O)(1\/sqrt(N))$. While this rate is independent of dimension, the variance $sigma^2$ often increases with dimension or payoff complexity, which can lead to low efficiency in practical settings.

@plain-monte-carlo-se-vs-sample-size illustrates this behaviour: the log-log plot of S.E. versus sample size $N$ follows a straight line with slope approximately $-0.5$, corroborating the theoretical rate.

#figure(
  caption: [Plain Monte Carlo: Standard Error vs Sample Size],
)[
  #image("assets/plain-monte-carlo-se-vs-sample-size.png",
    width: 10cm,)
] <plain-monte-carlo-se-vs-sample-size>

Achieving higher accuracy by increasing $N$ is computationally expensive: o reduce S.E. by a factor of tem requires roughly one hundred times more samples. This limitation motivates the development of _variance reduction techniques_, which aim to decrease the variance $sigma^2$ of the estimator without increasing the sample size.

#pagebreak()
= Variance Reduction Techniques
Although the plain Monte Carlo method provides unbiased estimates and is straightforward to implement, its slow convergence rate and potentially large variance often make it computationally inefficient for high-precision estimation. To address this, a range of variance reduction techniques have been developed to improve estimator accuracy without proportionally increasing the number of samples.

This section introduces several widely used methods: antithetic sampling, which exploits negative correlation between paired samples; control variates, which leverage known expectations of correlated quantities; control functionals, which use functional approximations within reproducing kernel Hilbert spaces; stratified sampling, which ensures balanced coverage of the sample space; and importance sampling, which concentrates samples in regions of high contribution to the integral.

Other advanced techniques, such as quasi-Monte Carlo methods and multi-level Monte Carlo, are often use to obtain further gains in efficiency but are beyond the scope of this report.

== Antithetic Sampling
// TODO: write this section
=== Examples
=== Limitations


== Control Variates
Suppose we can find another random variable $Y$ such that:

#enum(
  numbering: "i.",
  [$Y$ is correlated with $X$.],
  [Simulating $Y_i$ alongside $X_i$ is computationally inexpensive. ],
  [The expected value $EE[Y]$ is known exactly.]
)
Then we can construct a new estimator for $EE[X]$ that leverages the information from $Y$ to reduce variance.  The form of the new estimator is
$ hat(v)_"CV" = 1/N sum_(i = 1)^N X_i - b(1/N sum_(i = 1)^N Y_i - EE[Y]) $
and $Y$ is known as the _control variate_. Note that $hat(v)_"CV"$ is still an unbiased estimator of $EE[X]$ for any choice of $b in RR$.
Then the variance of $hat(H)$ is a quadratic in $b$:
$ "Var"[hat(H)] = 1 / N ["Var"[X] + b^2 "Var"[Y] - 2 b "Cov"(X, Y)], $
which is minimized when $b$ is chosen as
$ b^* = beta sigma_X / sigma_Y =  ("Cov"(X, Y)) / ("Var"(Y)) $
where $beta$ is the correlation between $X$ and $Y$. Intuitively, we are regressing $X$ against $Y$ and removing the component of $X$ that is explained by $Y$, resulting in lower variance#footnote[Usually, $b^*$ is estimated from a pilot run independent to the main simulation, otherwise bias may be introduced. In this view, $b^*$ is the estimated coefficient of $Y$ in the linear model $X ~ Y$.]. Under $b^*$, the variance of the new estimator is
$ "Var"[hat(v)_"CV"] = (1 - beta^2) "Var"[X] / n. $
This is lower than the variance of the plain Monte Carlo estimator by a factor of $1 - beta^2$. Therefore, the effectiveness of the control variate method depends on the strength of the correlation between $X$ and $Y$; stronger correlation leads to greater variance reduction.

=== Examples
#enum(
  [*European Call Option.* Suppose we want to estimate the price of a European call option with strike price $K$ and maturity $T$ using Monte Carlo simulation. The payoff of the option is given by $h(S_T) = (S_T - K)^+$, where $S_T$ is the stock price at maturity. We can use the stock price $S_T$ itself as a control variate since its expected value under the risk-neutral measure is known.
  
  $ EE[S_T] = S_0 e^(r T) $

  The control variate estimator is given by
  $ H = dash(X) - b(dash(Y) - EE[Y]) $
  where $X = e^(-r T) h(S_T)$ and $Y = S_T$. The optimal coefficient $b^*$ can be estimated from a pilot simulation.

  The basic Monte Carlo estimate is compared against the control variate estimate below.
  
  [table of results]
  ],
  [*Straddle Option.* A straddle option is a combination of a call and a put option with the same strike price and maturity. The payoff of a straddle option is
  $ h(S_T) = (S_T - K)^+ + (K - S_T)^+. $

  Similar to the case of a European call option, the terminal stock price $S_T$ can serve as a control variate to reduce the variance of Monte Carlo estimates. Consider the following parameter values:
  $ K = 50, r = 0.02, sigma = 0.2, " and" T = 1. $
  We summarize the results for $S_0 = 40, 50, 60$ and sample size of $N = 10,000$ in @straddle-option-results.

  #figure(
    caption: [Straddle Option: control variates versus plain Monte Carlo],
  )[
  #table(
    columns: 7,
    align: center,

    [Strike Price], table.cell(colspan: 2,align: center)[$S_0 = 40$], table.cell(colspan: 2,align: center)[$S_0 = 50$], table.cell(colspan: 2,align: center)[$S_0 = 60$],
    [], [Plain MC], [CV], [Plain MC], [CV], [Plain MC], [CV],
    [Estimate], [10.3312], [10.4381], [7.8795], [7.8737], [12.7905], [12.7598],
    [S.E.], [0.0614], [0.0400], [0.0628], [0.0593], [0.1018], [0.0416],
    [R.E.], [0.59%], [0.38%], [0.80%], [0.75%], [0.83%], [0.33%],
    [$hat(beta)^2$], table.cell(colspan: 2)[57.54%], table.cell(colspan: 2)[10.89%], table.cell(colspan: 2)[83.34%],
  )] <straddle-option-results>

  A variance reduction of approximately 57.5% is observed when using $S_T$ as a control variate for $S_0$. However, as $S_0$ approaches the strike price $K = 50$, the efficiency of this control variate diminishes, as illustrated in @vrf-S0.
  #figure(
    caption: [Variance reduction from the control variate as $S_0$ varies from 40 to 60, with $K = 50$.],
  )[
  #image("assets/control-variates-vrf-S0.png",
    width: 10cm,)
  ] <vrf-S0>

  When $S_0 = K = 50$, the variance reduction is approximately $10.9%$. This decline in performance can be explained by the weak correlation between $S_T$ and the straddle payoff $h(S_T)$. As shown in @straddle-option-linear-scatter, the relationship between these two quantities is highly non-linear when $S_0 approx K$. Because the payoff function is symmetric and exhibits a kink at $S_T = K$, the linear CV $S_T$ fails to capture much of the variation in $h(S_T)$, leading to suboptimal variance reduction.

  #figure(
    caption: [$h(S_T)$ vs. $S_T$ for straddle Option, $S_0 = K = 50$, using $S_T$ as a control variate.],
  )[
    #image("assets/straddle-option-linear-scatter.png",
      width: 10cm,)
  ] <straddle-option-linear-scatter>

  To improve upon this, additional control variates can be introduced. In particular, using both $S_T$ and $S_T^2$ may result in a better fit. Under the risk-neutral measure, the expected value of $S_T^2$ is known analytically as
  $ EE[S_T^2] = EE[S_T]^2 + "Var"(S_T) =  S_0^2 e^(2 r T + sigma^2 T). $
  The optimal coefficients $b_1^*$ and $b_2^*$ can then be estimated via multiple linear regression of $h(S_T)$ on $S_T$ and $S_T^2$. The new control variate estimator is given by
  $ hat(H) = dash(X) - b_1^*(dash(Y_1) - EE[Y_1]) - b_2^*(dash(Y_2) - EE[Y_2]) $
  where $Y_1 = S_T$ and $Y_2 = S_T^2$. 
  Incorporating both control variates yields a variance reduction of approximately 86.2%, a substantial improvement over the single control case.

  #figure(
    caption: [$h(S_T)$ vs. $S_T$ for straddle Option, $S_0 = K = 50$, using $S_T$ and $S_T^2$ as control variates.],
  )[
    #image("assets/straddle-option-quadratic-scatter.png",
      width: 10cm,)
  ] <straddle-option-quadratic-scatter>
  
  This result highlights that for payoffs exhibiting strong non-linearity such as the straddle option, augmenting the control variate basis with higher-order terms could significantly enhance variance reduction efficiency.],
)

=== Limitations
The requirement that $EE[Y]$ is known exactly is often prohibitive -- in many practical scenarios, it may not be feasible to find a control variate with a known expected value. There have been some attempts to get around this by designing a finite basis of control variates
$ Phi(Y) = {psi_1 (Y), psi_2 (Y), ..., psi_m (Y)} $
each with a known expectation $ EE[psi_i (Y)]$. This approach is powerful but inherently limited by the expressiveness of the finite basis. The space of functions that can be represented in this basis may not be rich enough to approximate the relationship between $X$ and $Y$ well. If the dependence is highly non-linear or complex, capturing it accurately would require a very large number of basis functions, making the method computationally infeasible.

To overcome these limitations, the control functional (CF) framework introduced in #cite(label("10.1111/rssb.12185")) places the payoff function within a _richer function space_, enabling flexible, non-parametric approximation of complex relationships.

== Control Functionals
The key idea of control functionals is to replace the _finite linear basis_ with a *function space* that can approximate a much broader class of relationships between random variables.

Suppose we wish to estimate $EE[h(Z)]$ where $Z$ is a random variable with known density $p(z)$ and $h$ is a given integrable function. In this setting, the classical control variate idea can be seen as constructing a function $g(Z)$ (the control variate) that approximates $h(Z)$ from basis functions $psi_i (Z)$. That is,
$ h(Z) approx g(Z) = b_1 psi_1 (Z) + b_2 psi_2 (Z) + ... + b_m psi_m (Z) $
where $EE[g(Z)]$ is known from $EE[psi_i (Z)]$ and the coefficients $b_i$ are chosen to minimize the mean squared error between $h(Z)$ and $g(Z)$ over the samples $Z_i$.

The control functional framework extends this idea by allowing $g$ to come from a _reproducing kernel hilbert space_ denoted $cal(H)_+$. This hilbert space is associated with a positive-definite kernel $k_+(dot.c, dot.c)$ that implicitly defines an infinite-dimensional set of basis functions $psi_i (z) = k_+ (z, z_i)$ for any $z_i in Omega$. So we can write
$ h(Z) approx g(Z) = sum_(i = 1)^(oo) b_i psi_i (Z) $
where $c in RR$ is a constant term. The kernel $k_+ (dot.c, dot.c)$ is constructed such that the corresponding functions have analytically known expectations (in particular, $EE[psi_i (Z)] = 1$ for all $i$).

Now the control functional $g in cal(H)_+$ is the function that best approximates $h$ in $cal(H)_+$. This can be realized as the solution to the following regularized least squares problem:
#math.equation(numbering: "(1)", block: true)[
$ g = arg min_(g in cal(H)_+) 1/N sum_(i = 1)^N (h(z_i) - g(z_i))^2 + lambda ||g||_(cal(H)_+)^2 $] <rls-eq>
where $lambda > 0$. After obtaining $g$, the control functional estimator for $EE[h(Z)]$ is given by
$ hat(H)_"cf" = dash(h) - dash(g) + EE[g (Z)] $
where $dash(h)$ and $dash(g)$ denote the sample means of $h(Z_i)$ and $g (Z_i)$ respectively, and $EE[g (Z)]$ is computed as the sum of $b_i$'s. Because $g$ resmembles $h$, the residuals $h(Z_i) - g(Z_i)$ have lower variance than $h(Z_i)$ alone, leading to a more efficient estimator.

=== Stein Operator
The construction of the desired RKHS $cal(H)_+$, we rely on the _Stein operator_ $cal(S)_p$, which encodes information about the
target distribution $p$. Fo a sufficiently smooth function $phi$, the Stein operator is defined as
$ cal(S)_p [phi](z) = nabla_z dot.c phi(z) + phi(z) dot.c nabla_z log p(z) $
where $nabla_z$ denotes the gradient with respect to $z$.

#enum(
  [*$cal(H) -> cal(H)_0.$* A key property of this operator is that, for any $phi$ belonging to an appropriate function class $cal(H)$,  
  $ EE[cal(S)_p [phi](Z)] = 0.$

  Under mild regularity conditions, we can define the hilbert space
  $ cal(H)_0 = { cal(S)_p [phi] : phi in cal(H)}, $
  where $cal(H)$ is a base RKHS with kernel $k(dot.c , dot.c)$. The induced kernel on $cal(H)_0$ is given by
  $ k_0 (z, z') = nabla_z dot.c nabla_(z') k(z, z') + u(x) dot.c nabla_(z') k(z, z') + u(z') dot.c nabla_z k(z, z') + u(z) dot.c u(z') k(z, z') $ 
  where $u(z) = nabla_z log p(z)$ is the _score function_ of the distribution.

  Equivalently, $k_0$ can be expressed more compactly as
  $ k_0 (z, z') &= cal(S)_p^z cal(S)_p^(z') [k(z, z')], $
  where $cal(S)_p^z$ and $cal(S)_p^(z')$ denote the Stein operator acting on the first and second arguments of $k$ respectively.
  This kernel satisfies the zero-mean property
  $ integral k_0 (z, z') p(z') dif z' = 0 $
  for almost all $z in Omega$,
  which ensures that any function $g in cal(H)_0$ has zero expectations under $p$, i.e., $EE_p [g(Z)] = 0$.
  This property is fundamental to control functionals as it guarantees that all functions 
  in $cal(H)_0$ can serve as valid, analytically tractable control variates.],
  [*$cal(H)_0 -> cal(H)_+. $* As established, the stein operator induces a RKHS $cal(H)_0$ of zero-mean functions.
  This is ideal for constructing mean-zero control variates, but to approximate a general target function $h$, we
  also need to represent its mean component. 
  
  To achieve this, we augment $cal(H)_0$ with the one-dimensional RKHS of constant functions:
  $ cal(C) = {c : c in RR}, $
  which has the reproducing kernel $k_C (z, z') = 1$ for all $z, z' in Omega$.
  The augmented RKHS is then defined as
  $ cal(H)_+ = cal(C) + cal(H)_0 = {c + g : c in RR, g in cal(H)_0}. $
  By the properties of RKHS addition, the corresponding kernel is simply the sum of the component kernels:
  $ k_+ (z, z') = k_C (z, z') + k_0 (z, z') = 1 + k_0 (z, z'). $
  This construction gives us a rich function space $cal(H)_+$ that can approximate a wide variety of target functions $h$, while still allowing for analytically tractable expectations via the Stein operator.
  ]
)

=== Example
Let $Z ~ N(0, 1)$, so $p(z) = (2 pi)^(-1/2) e^(-z^2\/2)$, and consider the standard Radial Basis Function (RBF) kernel:
$ k(z, z') = exp(-(z - z')^2\/2 alpha^2) $
that defines an RKHS $cal(H)$ of smooth functions, with basis functions ${k(z, z_i)}_(i=1)^n$ for any set of points ${z_i}_(i=1)^n$. The standard normal density has score function
$ u(z) = nabla_z log p(z) = -z,$
so the corresponding _Stein operator_ is
$ cal(S)_p [f] (z) = f'(z) - z f(z), $
where $f$ is differentiable. Applying the Stein operator to the RBF kernel, we obtain
$ k_0 (z, z') &= cal(S)_p^z cal(S)_p^(z') [k(z, z')] = e^(-(z - z')^2\/2 alpha^2) [1/alpha^2 - (z - z')^2 (1/alpha^4 + 1/alpha^2) + z z']. $
See @RKHS-kernel-functions for a visualization of elements from $cal(H)$ and $cal(H)_0$.

#figure(
  caption: [Elements of $cal(H)$ and $cal(H)_0$ for the standard normal distribution and RBF kernel with lengthscale $alpha = 1.0$],
)[
  #image("assets/RKHS-kernel-functions.png")
] <RKHS-kernel-functions>

#enum(
  [*Butterfly Multistrike Option.* Consider a deep out-of-the-money butterfly option with strike prices $K_1 = 90, K_2 = 100$, and $K_3 = 120$. The payoff of the option at maturity $T$ is given by
  $ h(S_T) = (S_T - K_1)^+ - 2 (S_T - K_2)^+ + (S_T - K_3)^+. $
  We compare the performance of the basic Monte Carlo estimator against the control functional estimator using the standard normal distribution and RBF kernel with $alpha = 1.0$. The parameters used are:
  $ r = 0.05, sigma = 0.3, T = 1. $
  For $S_0 = 60$ and a sample size of $N = 10,000$, we report the results in @butterfly-option-results below:

  #figure(
    caption: [Butterfly Option: control functionals versus control variates and plain Monte Carlo],
  )[
  #table(
    columns: 4,
    align: center,

    [], [Estimate], [S.E.], [R.E.],
    [Plain MC], [0.1585], [0.0179], [11.28%],
    [Control Variate], [0.1608], [0.0181], [11.27%],
    [Control Functional], [0.1628], [0.0023], [1.39%],
    [Theoretical], table.cell(colspan: 1)[0.1628]
  )
  ] <butterfly-option-results>

  The control functional estimator achieves a significant variance reduction compared to both the plain Monte Carlo and control variate methods. Specifically, the relative error is reduced from 11.28% to 1.39%. This demonstrates the effectiveness of control functionals in capturing complex relationships in the payoff function that linear control variates may miss.
  
  #figure(
    caption: [Butterfly Multistrike Option Pricing with $S_0 = 60$ and $N = 10,000$.],
  )[
    #image("assets/control-functionals-example.png",
      width: 16cm,)
  ] <control-functionals-example>
  ]
)

=== Limitations
While control functionals provide a powerful framework for variance reduction, they come with certain limitations. This section discusses several key limitations, focusing on computational cost, kernel selection, and theoretical assumptions.

#enum(
  numbering: "1.",
  [*Computational Cost.* A principal drawback of control functional methods lies in their computational complexity. Solving the RLS problem in equation (@rls-eq[]) involves inverting an $N times N$ kernel matrix $K$, which incurs a computational cost of $cal(O)(N^3)$ in time and $cal(O)(N^2)$ cost in memory.
  For large sample sizes $N$, the variance reduction benefits must be weighed against the  increased computational burden.
  In contexts where the cost of function evaluation is low, the overhead of control functionals may not be justified.],
  [*Kernel Selection.* The performance of control functionals is sensitive to the choice of kernel and its hyperparameters (e.g., lengthscale in RBF kernels). Selecting an appropriate kernel may require domain knowledge or empirical tuning, which can be non-trivial.],
  [*Theoretical Assumptions.* The theoretical appeal of control functionals stem from their super-root-$N$ rate of convergence, which, under ideal conditions, can outperform other variance reduction techniques. However, these rates are guaranteed only under specific regularity assumptions about the kernel and target function. In particular, this convergence property assumes that $ sup_(z in Omega) k_0(z, z) < oo, $ ensuring that the reproducing kernel $k_0$ is bounded. Commonly used kernels such as the RBF kernel
  $ k(z, z') = exp(-(z - z')^2\/2 alpha^2) $
  examined in the previous example, the derived kernel $k_0$ does not satisfy this assumption. To address this, alternative kernels with bounded $k_0$ can be employed. For instnace, Oates _et al._ #cite(label("10.1111/rssb.12185")) modifies the RBF kernel by introducing polynomial decay terms to ensure boundedness:
  $ k'(z, z') = (1 + alpha_1 (z^2 + z'^2))^(-1) exp(-(z - z')^2 \/ 2alpha_2^2) $
  where $alpha_1, alpha_2 > 0$.]
)

Despite these limitations, ongoing research continues to refine and extend the control functional framework, aiming to mitigate these challenges and broaden its applicability.


== Stratified Sampling
// TODO: write this section

==  Importance Sampling
Importance sampling (IS) reduces variance by sampling more frequently from regions of the domain that contribute most to the expectation. Instead of drawing samples from the original distribution $f(x)$, samples are taken from a proposal distribution $g(x)$ and reweighted to correct for the change in measure.

To estimate $v = EE[h(X)]$, the standard Monte Carlo estimator
$ hat(v) = 1/N sum_(i=1)^N h(X_i), quad X_i ~ f(x) $
is replaced by the importance sampling estimator
$ hat(v)_"is" = 1/N sum_(i=1)^N h(X_i) f(X_i) / g(X_i), quad X_i ~ g(x). $

The choice of $g(x)$ is critical: a good proposal distribution closely matches the shape of the function $|h(x)|f(x)$, leading to lower variance in the weighted estimator.

#pagebreak()
= Applications in Finance
// TODO: write this section

#pagebreak()
= Conclusion
// TODO: write this section

#pagebreak()
#bibliography("references.bib", full: true)

#pagebreak()
#include("appendix.typ")