// Capstone Report
// Author: Chin Zhe Ning
// Date: 2025-11-05

#import "utils.typ": *

// setup
#set page(margin: (x: 2.54cm, y: 2.54cm), numbering: "1")
#set heading(numbering: "1.")
#set par(first-line-indent: 1em, spacing: 1.5em,)
#set enum(indent: 1.25em, spacing: 1.25em)

#show heading: set text(
  size: 14pt,

  fill: blue.darken(70%))
#show heading: it => {
  v(0.5em)
  it
  v(0.5em)
}

#set text(size: 11pt, font: fonts.serif)
#set par(leading: 1.15em)  // line spacing

#set table(stroke: 0.5pt)

#show figure: it => {
  v(1em)
  it
  v(1em)
}

#set enum(indent: 1.5em, spacing: 1.5em)

#cover_page(
  title: "Monte Carlo Methods and their Applications to Finance",
  student: "Chin Zhe Ning",
  supervisor: "Dr. Ren Weiqing",
  examiner: "",
  module: "MA4198 Capstone Project",
  institution: "National University of Singapore",
  info: "Semester 1, AY2025/2026",
)

#pagebreak()
#outline(depth: 2)

#pagebreak()

= Introduction

This report is split into three main sections.
The first section lays out the basic concepts and definitions in option pricing theory, and how derivatives can be modeled as expectations via the no-arbitrage principle as presented by Wang (2012) in the book "Monte Carlo Simulation and Applications in Finance"  #cite(label("textbook")).
 We will introduce the plain Monte Carlo estimator for computing estimates of these expectations and also discuss its limitations.

In the second section, we first present several classical variance reduction techniques, including antithetic sampling and control variates, along with numerical examples to illustrate their effectiveness in option pricing problems. These methods impove the efficiency of Monte Carlo estimators by reducing their variance without increasing the number of samples. Most of this material is adapted from #cite(label("textbook")) and are well-established techniques in the field. We then introduces control functionals, a more recent advancement in variance reduction techniques based on reproducing kernel hilbert spaces and Stein's method as presented in #cite(label("cf")).

The last section focuses on the importance sampling technique, which is widely used in practice for pricing rare-event options. We discuss the theoretical foundations of importance sampling, including the choice of optimal proposal distributions to minimize estimator variance. Numerical examples are provided to demonstrate the effectiveness of importance sampling in reducing variance for pricing deep out-of-the-money options. Applications of importance sampling to modelling portfolio credit risk are also discussed.

== Option Pricing

Suppose that the payoff of a derivative is given by a function $h(Z)$ of a random
variable $Z$ related to the underlying asset. By the principle of no-arbitrage, the derivative's price can be modeled as the expected discounted payoff under the risk-neutral measure $QQ$
$ v = EE^QQ [e^(-r T) h(Z)], $

//$ v = EE[e^(-r T) h(X)], $
where $r$ is the _risk-free interest rate_ and $T$ is the time to maturity[cite], when the payoff is realized. There is a lot that can be said about the risk-neutral measure, but for our purposes, it suffices to know that it is a probability measure under which the discounted asset prices are _martingales_.
Under this measure, all assets are expected to grow at the risk-free rate, simplifying pricing to the computation of expected discounted payoffs.

Given an appropriate model for the underlying asset price dynamics, we can sample from the distribution of $Z$ and compute samples of $h(Z)$.
A basic model for the underlying asset price $S_t$ at time $t$, is geometric Brownian motion:
$ S_t = S_0 exp { (r - 1/2 sigma^2) t + sigma  W_t }, quad W_t ~ N(0, t) $
where $S_0$ is the intial asset price, $r$ is the risk-free interest rate, $sigma$ is the (constant) volatility, and $W_t$ is a standard Brownian motion.
Because of its simple theoretical properties, we will assume this model throughout the report.

=== Examples <option-examples>
#enum(
  numbering: "1.",
  [*European Call Option.* For an European call option with a strike price $K$ and underlying stock price $S_t$, the discounted payoff at time $T$ is
  $ h(S_T) = e^(-r T)(S_T - K)^+, $
  where $z^+ = max(z, 0)$. The price of the option is
  $ v = EE[e^(-r T ) (S_T - K)^+] $
  which can be solved analytically from the Black-Scholes formula #cite(label("black1973pricing")) as
  $ v = S_0 Phi(sigma sqrt(T) - theta) - K e^(-r T) Phi(-theta) $
  where the value of $theta$ is 
  $ theta = 1 / (sigma sqrt(T)) log K / S_0 + (sigma / 2 - r / sigma)sqrt(T). $],
  [*Asian Call Option.* For a discretely monitored average price call option, the
  payoff is given by
  $ h(S_(t_1), S_(t_2), ..., S_(t_m)) = e^(-r T )(1/m sum_(i= 1)^m S_(t_i) - K)^+, $
  where $0 < t_1 < dots.c < t_m = T$ are a fixed set of dates. Closed-form solutions for this path dependent option does not exist but for the case where the average is geometric, an analytical solution is available #cite(label("kemna1990pricing")). Hence, we resort to Monte Carlo methods to estimate the option price.],
)

These are two running examples that will be used throughout the report to illustrate the concepts and methods discussed.

== Plain Monte Carlo
Suppose that we are trying to estimate $v = EE[h(Z)]$, the option price (for ease of notation, we omit the discount factor $e^(-r T)$ in this section) where $X = h(Z)$ is a function of a random variable $Z$ with density $p_Z$. To generate samples $X_i = h(Z_i)$ we simulate $N$ independent and identically distributed (i.i.d.) samples ${Z_i}_(i=1)^N$.
The plain Monte Carlo estimate of $v$ is the sample average
$ hat(v) = 1/ N sum_(i=1)^(N) h(Z_i) = 1/ N sum_(i=1)^(N) X_i. $
This is an unbiased estimate $hat(v)$ of the true mean $v$, i.e., $EE[hat(v)] = v$,
with variance
$ "Var"(hat(v)) = 1/ N^2 sum_(i=1)^N "Var"(X_i)  = "Var"(X) / N, $
As the number of samples $N$ increases, the variance of $hat(v)$ decreases proportionally to $1\/N$. Consequently, the estimator 
$hat(v)$ converges to the true mean $v$ in the long run.

=== Examples <plain-monte-carlo-examples>
#enum(
  numbering: "1.",
  [*European Call Option.* Consider the problem of estimating the price of a European call option described in @option-examples with $K = 60$ and $T = 1$. Assuming the parameters
  $ S_0 = 50, r = 0.05, sigma = 0.2,$ we estimate the option price using the plain Monte Carlo method for various sample sizes $N$.
  The Monte Carlo estimate, standard error (S.E.), and relative error (R.E.) are presented in @plain-monte-carlo-results.

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
  [*Asian Call Option.* Consider the problem of estimating the price of a discretely monitored Asian call option described in @option-examples with $k = 60$ and $T = 1$. Assuming the parameters
  $ S_0 = 50, r = 0.05, sigma = 0.2 $
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

=== Limitations
The numerical results in the previous section for the European call and Asian options demonstrate that, the plain Monte Carlo estimator is consistent but it has a slow rate of convergence rate. This can be analyzed theoretically via the _Central Limit Theorem_.

Suppose that $"Var"(X) = sigma^2 < oo$, then
$ sqrt(N) (hat(v) - v) -> N(0, sigma^2), $
as $N -> oo$. This implies that  the standard error (S.E.) of $hat(v)$ decays at rate $cal(O)(N^(-1/2))$.

#figure(
  caption: [Plain Monte Carlo: standard error vs sample size.],
)[
  #image("assets/plain-monte-carlo-se-vs-sample-size.png",
    width: 12cm,)
] <plain-monte-carlo-se-vs-sample-size>

This is illustrated in @plain-monte-carlo-se-vs-sample-size where the log-log plot of S.E. versus sample size $N$ follows a line with slope approximately $-1\/2$ for both the European and Asian Call Option estimates, matching the theoretical rate. However, in many situations, this rate of convergence may be too slow to achieve the desired accuracy within a reasonable computational budget.

When the option becomes _deep out-of-the-money_, e.g., $K >> S_0$, the frequency of exercising is extremely low, so the distribution of the payoff $h(X)$ is very skewed. This leads to large $sigma^2$ and slow convergence. In @plain-monte-carlo-re-vs-strike-price, we vary the strike price of the European call option from $K = 40$ to $K = 120$ and observe the relative error (R.E.) on a fixed sample size of $N = 10,000$. As $K$ increases, the R.E. explodes resulting in worthless estimates.

#figure(
  caption: [Plain Monte Carlo: Relative Error vs Strike Price],
)[
  #image("assets/plain-monte-carlo-re-vs-strike-price.png",
    width: 12cm,)
] <plain-monte-carlo-re-vs-strike-price>

To get a more accurate estimate, increasing the sample size $N$ can become computationally prohibitive. Increasing the precision of the estimate by one decimal place (reducing S.E. by a factor of ten) requires one hundred times more samples. This motivates the usage of _variance reduction techniques_, to decrease the variance $sigma^2$ of the estimator without increasing the sample size.


#pagebreak()
= Variance Reduction Techniques

This section discusses two widely used methods in Monte Carlo to reduce estimator variance, antithetic sampling and control variates, following the presentation in #cite(label("textbook")).
Control functionals will also introduced as an extension to Control Variates as a more recent approach. 
All three techniques construct a new Monte Carlo estimator $hat(v)_"new"$ that is still unbiased but has lower variance than the plain estimator $hat(v)$ for the same sample size $N$. 


== Antithetic Sampling
Antithetic sampling is a varaince reduction technique that involves generating i.i.d. pairs of random variables $(X_i, Y_i)$ such that

#enum(
  numbering: "i.",
  [$X_i$ and $Y_i$ have the same distribution.],
  [$X_i$ and $Y_i$ are _negatively correlated_.]
)
We can construct a new unbiased estimator for $v = EE[X]$ of the form
$ hat(v)_"AS" = 1/N sum_(i = 1)^N (X_i + Y_i)/2. $
If $beta$ is the correlation between $X$ and $Y$, and $beta < 0$, then the variance of $hat(v)_"AS"$
$ "Var"[hat(v)_"AS"] = (1 + beta) "Var"[X] / (2N) $
is lower than the variance of the plain Monte Carlo estimator with $2N$ samples, that is,
$ "Var"[hat(v)] = "Var"[X]/(2N) $
Therefore, the variance reduction of the antithetic sampling method depends the strength of the negative correlation. If $beta << 0$, then variance reduction is good. In practice, the antithetic pairs $(X_i, Y_i)$ are often generated by inverting an underlying random variable. For example, if $X_i = h(U_i)$ where $U_i ~ "Uniform"(0, 1)$, then we can set $Y_i = h(1 - U_i)$. This construction ensures that $X_i$ and $Y_i$ are negatively correlated if $h$ is monotonic. To extend this to other random variables $Z$, we can compose $h$ with the _inverse cumulative distribution function_ of $Z$ as follows, $g = h compose F^(-1)_Z$, and generate the antithetic pair as
$ X_i = g (U_i), quad Y_i = g (1 - U_i). $





=== Examples
#enum(
  numbering: "1.",
  [*European Call Option.* The European call option described in @option-examples assumes  the underling asset price is a geometric Brownian motion, i.e.,  $S_T$ can be expressed as a function of $Z$:
  $ S_T = S_0 exp {(r - 1/2 sigma^2) T + sigma sqrt(T) Z}, quad Z ~N(0, 1) $
  This means that the payoff $h$ is a function of $Z ~ N(0, 1)$. Since $Z$ is symmetric around zero, the antithetic pair can be generated from $Z_i$ as $ X_i = h(Z_i), quad Y_i = h(-Z_i). $
  
  Simulating with the same parameters as in @plain-monte-carlo-examples and a sample size of $N = 10,000$, we compare the plain Monte Carlo estimator $hat(v)$ with the antithetic sampling estimator $hat(v)_"AS"$ in @european-call-antithetic-results for $S_0 = 40, 50, 60$.

  #figure(
    caption: [European Option: antithetic sampling versus plain Monte Carlo],
  )[
  #table(
    columns: 7,
    align: center,
    [Strike Price], table.cell(colspan: 2, align: center)[$S_0 = 40$], table.cell(colspan: 2,align: center)[$S_0 = 50$], table.cell(colspan: 2,align: center)[$S_0 = 60$],
    [], [Plain MC], [AS], [Plain MC], [AS], [Plain MC], [AS],
    [Estimate], [12.3593], [12.2832], [5.2996], [5.4026], [1.6731], [1.6009],
    [S.E.], [0.0950], [0.0233], [0.0731], [0.0385], [0.0442], [0.0285],
    [R.E.], [0.77%], [0.19%], [1.38%], [0.71%], [2.64%], [1.78%],
  )] <european-call-antithetic-results>
  The standard error of $hat(v)_"AS"$ is lower than standard error of $hat(v)$. 
  ],
  [*Straddle Option.* A straddle option is a combination of a call and a put option with the same strike price $K$ and maturity $T$. The payoff of a straddle option is
  $  h(S_T) = (S_T - K)^+ + (K - S_T)^+. $
  This payoff is not monotonic and so using the same antithetic sample construction as in the European call option may not yield negative correlation between the antithetic pairs.
  The simulation parameters are 
  $ K = 50, r = 0.02, sigma = 0.2, "and " T= 1. $
  and we report the results for $S_0 = 40, 50, 60$ with sample size of $N = 10,000$ in @straddle-option-antithetic-results.

  #figure(
    caption: "Straddle Option: antithetic sampling versus plain Monte Carlo"
  )[
  #table(
    columns: 7,
    align: center,

    [Strike Price], table.cell(colspan: 2, align: center)[$S_0 = 40$], table.cell(colspan: 2,align: center)[$S_0 = 50$], table.cell(colspan: 2,align: center)[$S_0 = 60$],
    [], [Plain MC], [AS], [Plain MC], [AS], [Plain MC], [AS],
    [Estimate], [10.3312], [10.4426], [7.9219], [8.2118], [12.8876], [12.7170],
    [S.E.], [0.0614], [0.0347], [0.0626], [0.0897], [0.1035], [0.0702],
    [R.E.], [0.59%], [0.33%], [0.79%], [1.09%], [0.80%], [0.55%],
  )] <straddle-option-antithetic-results>

  For $S_0 = 40$, the relative error of $hat(v)_"AS"$ decreases ($0.24%$) is lower than the relative error of $hat(v)$, indicating that antithetic sampling is effective in reducing variance. However, as $S_0$ approaches the $K$, the relative error of $hat(v)_("AS")$ increases until it is higher than the relative error of $hat(v)$ as shown in @straddle-option-antithetic-re-vs-S_0.

  #figure(
    caption: [Straddle Option: Relative error of $hat(v)$ and $hat(v)_"AS"$ for $S_0 in [40, 60]$ ]
  )[
    #image(
      "assets/straddle-option-re-vs-initial-stock-price.png",
      width: 14cm
    )
  ] <straddle-option-antithetic-re-vs-S_0>
  ]
)

=== Limitations
The effectiveness of antithetic sampling relies on inducing a strong negative correlation $beta << 0$ between paired samples $(X_i, Y_i)$.
This technique works well when the payoff function 
$h$ is monotonic in the underlying random variable (e.g., European call option). In such cases, the antithetic sample $h(-Z)$ tends to move inversely to 
$h(Z)$, producing $beta << 0$.
However, for non-monotonic payoff functions, e.g. straddle option, it is possible for $beta > 0$.

When $S_0 -> K$, the values of $X_i$ and $Y_i$ tend to be similar, producing strong positive correlation $beta >> 0$. This is shown in @straddle-option-antithetic-scatter.


#figure(
  caption: [Scatter plot of antithetic pairs $Y_i$ vs. $X_i$ for straddle option, with $S_0 = K = 50$.],
)[
  #image("assets/straddle-option-antithetic-scatter.png",
    width: 12cm,)
] <straddle-option-antithetic-scatter>



== Control Variates
Suppose we can find another random variable $Y$ such that:

#enum(
  numbering: "i.",
  [$Y$ is correlated with $X$.],
  [Simulating $Y_i$ alongside $X_i$ is computationally inexpensive. ],
  [The expected value $EE[Y]$ is known exactly.]
)
We can construct a new estimator $hat(v)_"CV"$ for $v = EE[X]$ that leverages the information from $Y$ about $X$ to reduce the variance. The form of the new estimator is
$ hat(v)_"CV" = dash(X) - b(dash(Y) - EE[Y]), $
where $b in RR$ is a constant to be determined and $dash(X)$, $dash(Y)$ are the sample means of $X$ and $Y$ respectively. The variable $X$ is the _target variable_ and $Y$ is known as the _control variate_.
For any choice of $b in RR$, $hat(v)_"CV"$ is an unbiased estimate and the variance is a quadratic in $b$:
$ "Var"[hat(v)_"CV"] = 1 / N ("Var"[X] + b^2 "Var"[Y] - 2 b "Cov"[X, Y]), $
which is minimized when $b$ is chosen as
$ b^* = beta sigma_X / sigma_Y =  ("Cov"[X, Y]) / ("Var"[Y]) $
where $beta$ is the correlation between $X$ and $Y$. Intuitively, we are regressing $X$ against $Y$ and removing the component of $X$ that is explained by $Y$, resulting in lower variance#footnote[Usually, $b^*$ is estimated from a pilot run independent to the main simulation, otherwise bias may be introduced. In this view, $b^*$ is the estimated coefficient of $Y$ in the linear model $X ~ Y$.]. Under $b^*$, the variance of the new estimator is
$ "Var"[hat(v)_"CV"] = (1 - beta^2) "Var"[X] / n = (1 - beta^2) "Var"[hat(v)], $
which is lower than the variance of the plain Monte Carlo estimator by a factor of $1 - beta^2$. Therefore, the effectiveness of the control variate depends on $0 <= |beta| <= 1$.

=== Examples <control-variates-examples>
#enum(
  [*Asian Option.* Recall that the discounted payoff of our discretely monitored Asian call option is
  $ X = h(S_t_1, S_t_2, ..., S_t_m) = e^(-r T)(dash(S) - K)^+, quad dash(S) = 1/m sum_(i = 1)^m S_(t_i). $
  We can use the discretely monitored geometric average price call option as a control variate because its mean is known analytically #cite(label("kemna1990pricing")). The control variate has the form
  $ Y =  e^(-r T)(dash(S)_G - K)^+, quad dash(S)_G = (product_(i = 1)^m S_(t_i))^(1\/m) $
  with mean $EE[Y]$ has the analytical form
  $ dash(mu) = log S_0 + (r - sigma^2 / 2) 1/m sum_(i = 1)^(m)  t_i, quad  dash(sigma)^2 = sigma^2 / m^2 sum_(i = 1)^(m) (2m - 2i + 1) t_i. $
  $ EE[Y] = e^(-r T) [e^(dash(mu) + 1/2 dash(sigma)^2) Phi(dash(sigma) - theta) + K Phi(-theta)], quad theta = (log K - dash(mu))/dash(sigma). $

  Simulating with parameters,
  $ S_0 = 50, r = 0.05, sigma = 0.2, T = 1 $
  with sample size $N = 10,000$ and monitoring at $m = 100$ equally space time points, we estimate the Asian option price using both the plain Monte Carlo estimator $hat(v)$ and the control variate estimator $hat(v)_"CV"$ for $K = 60, 70, 80$ in @asian-option-control-variate-results.
  The plain Monte Carlo estimate $hat(v)$ is compared against the control variate estimate $hat(v)_"CV"$ below.
  
  #figure(
    caption: [Asian Option: Control variate versus plain Monte Carlo],
  )[
    #table(
      columns: 7,
    align: center,

    [Strike Price], table.cell(colspan: 2, align: center)[$S_0 = 60$], table.cell(colspan: 2,align: center)[$S_0 = 70$], table.cell(colspan: 2,align: center)[$S_0 = 80$],
    [], [Plain MC], [CV], [Plain MC], [CV], [Plain MC], [CV],
    [Estimate], [0.2887], [0.2720], [0.0159], [0.0107], [0.0000], [0.0009],
    [S.E.], [0.0130], [0.0009], [0.0032], [0.0005], [0.0000], [0.0001],
    [R.E.], [4.51%], [0.33%], [19.83%], [4.26%], [nan%], [13.39%],
    [$hat(beta^2)$], table.cell(colspan: 2)[99.53%], table.cell(colspan: 2)[97.92%], table.cell(colspan: 2)[-],
    )
  ] <asian-option-control-variate-results>

  As $K$ increases, the asian option becomes deep out-of-the-money, resulting in high relative error (R.E.) for the plain Monte Carlo estimator $hat(v)$. However, using the geometric average price call option as a control variate significantly reduces the variance of the estimator.
  For $K = 80$, $hat(v)$ fails to provide a meaningful estimate while $hat(v)_"CV"$ still produces a reasonable estimate with R.E. of $13.39%$. The high value of $hat(beta)^2$ indicates a strong correlation between the geometric and arithmetic average price call option payoffs, leading to effective variance reduction.

  ],
  [*Straddle Option.* Again, we use the terminal stock price $S_T$ to serve as a control variate. With the same parameters as in the antithetic sampling example, we simulate the straddle option and compare the plain Monte Carlo estimator with the control variate estimator in @straddle-option-results.

  #figure(
    caption: [Straddle Option with the linear control variate $Y = S_T$.],
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

  A variance reduction of approximately 57.5% is observed when using $S_T$ as a control variate for $S_0$. However, as $S_0$ approaches the strike price $K = 50$, the variance reduction for this control variate diminishes to 10.9%.

 This decline in performance can be explained by the weak correlation between $S_T$ and the straddle payoff $h(S_T)$. As shown in @straddle-option-linear-scatter, the relationship between these two quantities is highly non-linear when $S_0 approx K$. Therefore, the linear control variate $Y = S_T$ fails to capture much of the variation in $h(S_T)$, leading to suboptimal variance reduction. This shows that the same control variate may not be effective across different parameter regimes.

  #figure(
    caption: [$h(S_T)$ vs. $S_T$ for straddle Option, $S_0 = K = 50$, using $S_T$ as a control variate.],
  )[
    #image("assets/straddle-option-linear-scatter.png",
      width: 10cm,)
  ] <straddle-option-linear-scatter>

  To improve upon this, additional control variates can be introduced. By observation, the $h(S_T)$ could be better fitted by a quadratic function of $S_T$.
  Since the expected value of $S_T^2$ is known analytically
  $ EE[S_T^2] = EE[S_T]^2 + "Var"(S_T) =  S_0^2 e^(2 r T + sigma^2 T), $
  we can use it as a second control variate.
  The optimal coefficients $b_1^*$ and $b_2^*$ can then be estimated via the linear regression of $h(S_T) ~ S_T + S_T^2$. The new control variate estimator is given by
  $ hat(v)_"CV" = dash(X) - b_1^*(dash(Y_1) - EE[Y_1]) - b_2^*(dash(Y_2) - EE[Y_2]) $
  where $Y_1 = S_T$ and $Y_2 = S_T^2$. 
  Incorporating both control variates yields a variance reduction of approximately 86.2%, a substantial improvement over the linear control variate.

  #figure(
    caption: [Straddle Option: $h(S_T)$ versus $S_T$ with $S_0 = K = 50$, using $S_T$ and $S_T^2$ as control variates.],
  )[
    #image("assets/straddle-option-quadratic-scatter.png",
      width: 10cm,)
  ] <straddle-option-quadratic-scatter>
  
  This result highlights that augmenting the control variate with higher-order terms could significantly enhance variance reduction efficiency.],
)

=== Limitations
The requirement that $EE[Y]$ is known exactly is often prohibitive. In many practical scenarios, finding a control variate with a known expected value is hard. Even then, $Y$ may not be well correlated with $X$. To address this, some approaches introduce a feature map $Phi : RR -> RR^m$,
$ Phi(y) = (psi_1 (y), psi_2 (y), ..., psi_m (y))^T $
where each $psi_i$ has a known expectation under the distribution of $Y$, i.e. $ EE[psi_i (Y)]$.
This allows one to construct control variates as _linear combinations_ of multiple known features
$ Y' = sum_(i = 1)^m b_i  psi_i (Y) $
where $b_i$ are coefficients chosen to minimize the variance of the new estimator
$ hat(v)_"CV" = dash(X) - (dash(Y)' - EE[Y']). $
Still, the choice of the functions $psi_i$ is nontrivial as the set of functions that can be represented as linear combinations of $psi_i$ may not be rich enough to capture the true relationship between $X$ and $Y$ well.


== Control Functionals
The key idea of the control functional framework is to replace the set of linear combinations of $psi_i$ with a _rich function space_ that houses a wider variety of functions.
 
Suppose we wish to estimate $v = EE[h(Z)]$ where $Z$ is a random variable with known density $p(z)$ and $h$ is a given integrable function. 
In this setting, classical control variates construct a function $g(z)$ that approximates $h(z)$ as a linear combination of _basis functions_ $psi_i (z)$. That is,
$ h(z) approx g(z) = b_1 psi_1 (z) + b_2 psi_2 (z) + dots.c + b_m psi_m (z) $
where $EE[g(Z)]$ is computed from $EE[psi_i (Z)]$ and coefficients $b_i$ are chosen to minimize the _empirical risk_ between $h(z)$ and $g(z)$ over the samples $Z_i$. This approach is equivalent to performing a linear regression of $h(Z)$ against the basis functions $psi_i (Z)$ and minimizing the variance of $hat(v)_"CV"$.

Control functionals approximates $h(z)$ by allowing $g$ to reside in a _reproducing kernel hilbert space_ (RKHS), denoted $cal(H)_+$, associated with a positive-definite kernel $k_+(dot.c, dot.c)$
$ h(z) approx g(z) in cal(H)_+ >> "span"(psi_i, ..., psi_m). $
The RKHS is much richer than the linear span of finite basis functions, allowing $g$ to finely approximate $h$.
To obtain $g in cal(H)_+$, we minimize the _empirical risk_ in the regularized least squares problem:
#math.equation(numbering: "(1)", block: true)[
$g^* =  arg min_(g in cal(H)_+) 1/N sum_(i = 1)^N ||h(z_i) - g(z_i)||^2 + lambda ||g||_(cal(H)_+)^2 $] <rls-eq>
where $lambda > 0$. From the _Representer Theorem_, the optimal solution $g^*$ to (@rls-eq[]) admits a representation of the form 
#math.equation(numbering: "(1)", block: true)[
$ g^*(z) = sum_(i = 1)^N b_i k_+(z, Z_i) $] <rep-thm-eq>
and so (@rls-eq[]) reduces to solving a linear system for $b_i$ #cite(label("ScholkopfSmola2002")).
After obtaining $g^*$, the control functional estimato is given by
$ hat(v)_"CF" = 1/N sum_(i = 1)^N h(Z_i) - g^*(Z_i) + EE[g^*(Z)]. $
To ensure that $EE[g^* (Z)]$ can be computed from (@rep-thm-eq[]), the value $EE[k_+(Z, Z_i)]$ must be known for any fixed $Z_i$. This is achieved using the _Stein operator_, which can construct $k_+(dot.c, dot.c)$ from a base kernel $k(dot.c, dot.c)$.

=== Stein Operator
Given a sufficiently smooth density $p : Omega -> RR_(>=0)$, the Stein operator $cal(S)_p$ is defined as
$ cal(S)_p [phi](z) = nabla_z dot.c phi(z) + phi(z) dot.c u(z) $
where $nabla_z$ denotes the gradient with respect to $z$ and $u(z) = nabla_z log p(z)$ is the _score function_ of the distribution. _Stein's identity_ states that for any function $phi$ satisfying appropriate boundary conditions#cite(label("cf")), we have the zero-mean property
$ EE[cal(S)_p [phi](Z)] = 0. $
This allows us to construct a reproducing kernel $k_0 (dot.c, dot.c)$, by applying $cal(S)_p$ to both arguments of a base reproducing kernel $k(dot.c, dot.c)$, with the zero-mean property, that is,
$ EE[k_0 (Z, Z_i)] = 0, $
for almost all $Z_i in Omega$. The details of this construction and further constructing $k_+$ from $k_0$ are not relevant to the report and thus omitted for clarity. The interested reader is referred to #cite(label("cf")) for more details. See @RKHS-kernel-functions for examples representative elements from $cal(H)$ and corresponding elements from $cal(H)_0$ obtained via the Stein operator for the standard normal distribution and a modified Radial Basis Function (RBF) kernel.
$ k(z, z') = (1 + alpha_1 ||z||^2 + alpha_1 ||z'||^2))^(-1) exp(-||z - z||^2\/2 alpha_2^2) $
where $alpha_1, alpha_2 > 0$ are the decay scale and length scale terms respectively.

#figure(
  caption: [Elements of $cal(H)$ and $cal(H)_0$ for the standard normal distribution and modified RBF kernel with decay scale $alpha_1 = 0.001$ and length scale $alpha_2= 1.0$.],
)[
  #image(
    "assets/RKHS-kernel-functions.png",
    width: 20cm,
  )
] <RKHS-kernel-functions>

=== Methodology <cf-methodology>
The samples ${Z_i}_(i=1)^N$ drawn i.i.d. from $p(z)$ are split into two disjoint sets: the _training set_ $cal(D)_0 = {Z_i}_(i=1)^M$ of size $M$ and the _test set_ $cal(D)_1 = {Z_i}_(i=M+1)^N$. The training set is used to learn the control functional $g in cal(H)_+$ by solving the RLS problem in (@rls-eq[]). And then both sets are used to compute the control functional estimator $hat(v)_"CF"$. This split-sample approach is necessary to ensure that the estimator $hat(v)_"CF"$ is unbiased#footnote[If the same samples are used to both fit $g$ and compute $hat(v)_"CF"$, then bias may be introduced because $g$ is fitted to minimize the empirical risk on those samples.]. The form of the control functional estimator is
$ hat(v)_"CF" =  1/(N - M) bold(1)^T (bold("h")_1 - hat(bold("h"))_1) + (bold(1)^T (bold("K")_0 + lambda M bold(I))^(-1) bold("h")_0) / (1 + bold(1)^T (bold("K")_0 + lambda M bold(I))^(-1) bold(1)) $
where $bold("h")_0 = [h(Z_1), ..., h(Z_M)]^T$, $bold("h")_1 = [h(Z_(M+1)), ..., h(Z_N)]^T$ are the payoffs on the training and test set respectively. And  $bold(1) = [1, ..., 1]^T$, and $(bold("K")_0)_(i, j) = k_0(Z_i, Z_j)$ is the kernel matrix on the training set. The predictions for the test set payoffs $bold("h")_1$ is 
$ hat(bold("h"))_1 = bold("K")_(1, 0) (bold("K")_0+ lambda M bold(I))^(-1) bold("h")_0 + (bold(1) - bold("K")_(1, 0)(bold("K")_0 + lambda M bold(I))^(-1) bold(1))((bold(1)^T (bold("K")_0 + lambda M bold(I))^(-1) bold("h")_0) / (1 + bold(1)^T (bold("K")_0 + lambda M bold(I))^(-1) bold(1))), $
with $(bold("K")_(1, 0))_(i, j) = k_0(Z_(M + i), Z_j)$.

The size of the training set $M$ is chosen proportional to $N$ such that $M/N -> eta in (0, 1)$ and the regularization parameter $lambda$ is chosen such that $lambda = O(M^(-1/2))$. Under some regularity conditions on the kernel $k_0(dot.c, dot.c)$, the density $p$, it can be shown that the control functional estimator achieves _super-root-$N$ convergence_. Intuitively, this is because the control functional $g$ approximates $h$ more finely as $N$ (and hence $M$) increases. The more samples used to fit $g$, the lower the variance of the residual $h(Z) - g(Z)$, leading to faster convergence of $hat(v)_"CF"$.

=== Examples <control-functionals-example>
*Straddle Option.* We revisit the same straddle option example from @control-variates-examples, using the same parameters and sample size of $N = 10,000$. The performance of the plain Monte Carlo estimator is compared against the control functional estimator using the standard normal distribution and a modified RBF kernel with $alpha_1 = 0.001, alpha_2 = 0.25$, and $eta = 0.4$. The results for $S_0 = 50$ are summarized in @straddle-option-cf-results below:

#figure(
  caption: [Straddle Option: control functionals versus control variates and plain Monte Carlo],
)[
  #table(
    columns: 4,
    align: center,

    [], [Estimate], [S.E.], [R.E.],
    [Plain MC], [7.8795], [0.0628], [0.80%],
    [Control Variate], [7.8737], [0.0593], [0.75%],
    [Control Functional], [7.9251], [0.0055], [0.07%],
    [Theoretical], table.cell(colspan: 1)[7.9260]
  )
] <straddle-option-cf-results>

The control functional fits well to the non-linear payoff function of the straddle option as shown in @straddle-option-control-functionals-example, and the estimator achieves a significant variance reduction compared to both the plain Monte Carlo and control variate methods. Specifically, the relative error is reduced from 0.80% to 0.07%. This demonstrates the effectiveness of control functionals in capturing complex relationships in the payoff function that linear control variates may miss.

For deep-out-of-the-money options, control functionals can be combined with importance sampling to further enhance variance reduction, which will be discussed in the next section.

#figure(
  caption: [Straddle Option Pricing with $S_0 = 50$ and $N = 10,000$, using control functionals.],
)[
  #image("assets/straddle-option-control-functionals-example.png",
    width: 20cm,)
] <straddle-option-control-functionals-example>

It was claimed that under certain conditions, control functionals can achieve super-root-$N$ convergence. To empirically verify this, we simulate the straddle option price with varying sample sizes $N$ and compute the standard error of both the plain Monte Carlo estimator and the control functional estimator. The results are shown in @straddle-option-cf-convergence.

#figure(
  caption: [Straddle Option: Convergence of control functional versus control variate and plain Monte Carlo estimators],
)[
  #image("assets/control-functional-se-vs-sample-size.png",
    width: 14cm,)
] <straddle-option-cf-convergence>

The standard error of the plain Monte Carlo estimator decreases at the expected rate of $O(N^(-1\/2))$, as indicated by the slope of approximately -0.5 on the log-log plot. In contrast, the control functional estimator exhibits a steeper decline in standard error, with an estimated slope of approximately -0.9. This demonstrates that the control functional estimator achieves a convergence rate faster than $O(N^(-1\/2))$.


=== Limitations
Control functionals have beautiful theoretical properties and are a first step towards constructing more general non-parametric control variates, they are not without drawbacks. This section outlines some of these drawbacks, focusing on computational cost, kernel selection, and theoretical assumptions.

#enum(
  numbering: "1.",
  [*Computational Cost.* Solving the RLS problem in (@rls-eq[]) according to the steps in @cf-methodology requires inverting an $M times M$ kernel matrix $K$, which incurs a computational cost of $cal(O)(M^3)$ in time and $cal(O)(M^2)$ cost in memory. This scales poorly and is entirely infeasible for $M >> 1e^5$.
  For training sample sizes $M$, the variance reduction benefits must be weighed against the increased computational burden.
  In contexts where the cost of function evaluation is low, the overhead of control functionals may not be justified.],
  [*Kernel Selection.* The performance of control functionals is sensitive to the choice of kernel and its hyperparameters (e.g., $alpha_1, alpha_2$ in the modified RBF kernel) and the regularization parameter $lambda$.
  Selecting an appropriate kernel requires empirical tuning, incurring additional compute.],
  [*Theoretical Assumptions.* The theoretical guarantee of super-root-$N$ rate of convergence hinges on several assumptions about the smoothness of $h$, the density $p$, and the properties of the RKHS induced by the kernel $k_0$. 
  One assumption is that the reproducing kernel $k_0$ is bounded $ sup_(z in Omega) k_0(z, z) < oo. $ And this may not hold for common kernels, e.g., standard RBF kernel
  $ k(z, z') = exp(-||z - z'||^2\/2 alpha^2). $
  To address this, alternative kernels with bounded $k_0$ can be employed. For instance, Oates et al. #cite(label("cf")) proposes the modified RBF kernel by introducing polynomial decay terms to ensure boundedness:
  $ k(z, z') = (1 + alpha_1 ||z||^2 + alpha_1 ||z'||^2)^(-1) exp(-||z - z'||^2 \/ 2alpha_2^2) $
  where $alpha_1, alpha_2 > 0$. The $k_0$ induced by this modified kernel is bounded as shown in Appendix A.
  ]
)

Despite these limitations, ongoing research continues to refine and extend the control functional approach, aiming to mitigate these challenges and broaden its applicability. One prominent direction replaces the fixed RKHS representation of control functionals with neural network–based models, leveraging the flexibility of deep architectures to handle high-dimensional and non-smooth integrands. For instance, #cite(label("muller2020neural")) introduced Neural Control Variates, a generalization of the control functional approach that parameterizes the control variate using normalizing flows, thereby enabling expressive, invertible transformations of the sampling distribution while preserving unbiasedness through Stein-type constraints. Building on this foundation, other works #cite(label("si2020scalable")), #cite(label("sun2023meta")) proposed scalable stochastic optimization and meta-learning formulations that unify kernel-, polynomial-, and neural-based control variates under a common variational or Stein-operator framework. 

#pagebreak()
=  Importance Sampling
Importance sampling (IS) reduces variance by sampling more frequently from regions of the domain that contribute most to the expectation.
== Introduction
Instead of drawing samples from the original distribution $f$, samples are taken from an alternative distribution $g$ and reweighted to correct for the change in measure.
The standard Monte Carlo estimator
$ hat(v) = 1/N sum_(i=1)^N h(X_i), quad X_i ~ f(x) $
is replaced by the importance sampling estimator
$ hat(v)_"IS" = 1/N sum_(i=1)^N h(X_i) f(X_i) / g(X_i), quad X_i ~ g(x) $
where the quantity $f(X_i)\/g(X_i)$ is known as the likelihood ratio or importance weight. 

Importance sampling has found widespread applications in finance, particularly in estimating rare-event probabilities such as Value-at-Risk (VaR) and Expected Shortfall (ES) where accurate estimation is difficult using standard Monte Carlo due to the low frequency of extreme losses. See @applications-in-finance for an illustration.

=== Variance Reduction
The main challenge in IS is to select an alternative distribution $g$ that reduces the variance of the estimator. Intuitively, the variance is minimized when more samples are drawn from regions where $h(x) f(x)$ is large.

#figure(
  caption: [Original distribution $f$ and two alternative distributions $g_1$ and $g_2$ for importance sampling.],
)[
]

The theoretically optimal choice $g^*(x)$ is proportional to $|h(x)|f(x)$, i.e.,
$ g^*(x) = c h(x) f(x), quad c in RR quad ==> quad 1/c = integral_(RR) h(x) f(x) dif x = v $
Under this ideal distribution, the variance of the importance sampling estimator is
$ "Var"[hat(v)_"IS"] = 1/N (integral_(RR) h(x)^2 (f(x)^2)/(g^*(x))dif x  - v^2) = 1/N (1/c integral_(RR) h(x) f(x) dif x - v^2) = 0. $
In practice, this is not feasible since it requires knowledge of $c$, which depends on the true mean $v$. But this ideal case provides useful insight into how to choose $g$, namely, $g$ should mimic the shape of $|h(x)| f(x)$ to reduce variance.

It is convenient to restrict $g$ to a parametric family of distributions
$ cal(G) = {g (x; theta) : theta in RR^d} $ that are easy to sample from. The goal then is to find an optimal $theta^*$ such that $g(x; theta^*)$ best approximates the ideal zero-variance distribution. This can be formalized as minimizing the Kullback-Leibler (KL) divergence between the ideal zero-variance distribution $g^*$ and the candidate $g(x; theta)$.
$ min_(g in cal(G)) integral_RR log (g^*(x)) / (g(x; theta)) dot.c g^*(x) dif x $
which is equivalent to the maximization problem
#math.equation(numbering: "(1)", block: true)[
$ max_(g in cal(G)) integral_RR  h(x) f(x) dot.c log g(x; theta) dif x. $] <kl-loss-eq>

This can be solved (approximately) in various ways depending on the parametric family $cal(G)$ and payoff function.

=== Mode Matching
A practical heuristic to choose $g(x; theta)$ is mode matching, which aligns the mode of $g$ with the mode of the target distribution $|h(x)| f(x)$.
Following the intuition from the ideal zero-variance distribution $g^*$, mass of $g$ is shifted towards the region where $h(x) f(x)$ is large, encouraging more frequent sampling from these high contribution regions.

For example, if we are choosing $g$ from the family of normal distributions $ cal(G) = {N(mu, 1) : mu in RR}, $ we can set $mu$ to be the mode of $|h(x)| f(x)$: 
$ mu^* = arg max_(x in RR) |h(x)|f(x). $
This simple approach can be effective for unimodal target distributions $|h(x)| f(x)$ in choosing an _approximately optimal_ solution to the maximization problem in (@kl-loss-eq[]). See @mode-matching-illustration for an illustration of mode matching in estimating the rare-event probability $P(Z >= b)$ for $b = 2.0$ where $Z ~ N(0, 1)$. The alternative distribution $g$ is obtained by shifting to $mu^* = 2.0$.


#figure(
  caption: [Under $N(0, 1)$, few samples fall in region $A = {Z >= b}$ and contribute to the estimate. Under $N(2, 1)$, more samples fall in $A$, and $g$ mimics $h(x)f(x)$ better.],
)[
  #image(
    "assets/importance-sampling-illustration.png",
    width: 16cm,
  )
] <mode-matching-illustration>

=== Iterative Cross Entropy Method
Instead of heuristically choosing $theta$ via mode matching, we can maximize the objective in (@kl-loss-eq[])  by iteratively updating $theta$. The gradient of the objective function with respect to $theta$ is given by
$ (partial)/(partial theta) L(theta) = integral_RR h(x) f(x) (partial )/(partial theta) log g(x; theta) dif x. $
This integral can be approximated using Monte Carlo by drawing samples $X_i ~ g(x; theta)$:
$ (partial)/(partial theta) L(theta) approx 1/N sum_(i=1)^N h(X_i) f(X_i) / g(X_i; theta) dot.c (partial)/(partial theta) log g(X_i; theta). $
If $L(theta)$ is maximized at $theta^*$ and is differentiable, with one global maximum, then the gradient at $theta^*$ is equal to zero. This suggests an iterative scheme where at each iteration $k$, we update $hat(theta)_k -> hat(theta)_(k+1)$ by solving for $theta$ in the equation
$ 1/N sum_(i=1)^N h(X_i) ell_(hat(theta)_k)(X_i) dot.c (partial)/(partial theta) log g(X_i; theta) = 0 $
where $ell_(hat(theta)_k)(X_i) = f(X_i) \/ g(X_i; hat(theta)_k)$ is the likelihood ratio evaluated at the current parameter estimate $hat(theta)_k$. 
In certain cases, this equation can be solved in closed form. For example, if $g$ is chosen from the family of normal distributions $cal(G) = {N(mu, I_m) : mu in RR^m}$, the update for $mu$ is simply
$ hat(mu)_(k+1) = (sum_(i=1)^N h(X_i) ell_(hat(mu)_k)(X_i) dot.c X_i) / (sum_(i=1)^N h(X_i) ell_(hat(mu)_k)(X_i) ) = (sum_(i=1)^N h(X_i) e^(-angle.l hat(mu)_k, X_i angle.r) dot.c X_i) / (sum_(i=1)^N h(X_i) e^(-angle.l hat(mu)_k, X_i angle.r) ). $
This is known as the general iterative cross-entropy (CE) method.

For more complicated parametric families $cal(G)$ where closed-form updates are not available, other optimization techniques such as stochastic gradient ascent, or EM-algorithm can be employed to maximize the objective in (@kl-loss-eq[]).

=== Examples
Similar to the problem estimating rare-event probabilities, deep out-of-the-money options are challenging to price accurately using standard Monte Carlo due to the low probability of exercise.
Importance sampling can be employed to focus sampling efforts on the tail regions of the underlying asset's distribution where the option is more likely to be in-the-money at maturity.
#enum(
  [*European Option.* Revisiting the European call option, but this time we consider a deep out-of-the-money case with parameters:
  $ S_0 = 50, r = 0.05, sigma = 0.2, T = 1. $
  We compare the performance of the plain Monte Carlo estimator against an importance sampling estimator across $K = 80, 100, 120$
  The alternative distribution $g$ is chosen via mode matching from the family of normal distributions $cal(G) = {N(mu, 1) : mu in RR}$.
  Taking the derivative of $h(x) f(x)$ and setting it to zero, we find that $mu^*$ satisfies
  $ S_0 exp {-1/2 sigma^2 + sigma sqrt(T) mu^*} (sigma sqrt(T) - mu^*) + e^(-r T) K mu^* = 0. $
  this is solved numerically using the bisection method yielding $mu^* approx 4.47$. We compare the result of the plain Monte Carlo estimate with the importance sampling estimate under $N(mu^*, 1)$ in @european-call-option-results.

  #figure(
    caption: [European Call Option: importance sampling versus plain Monte Carlo.],
  )[
    #table(
    columns: 7,
    align: center,

    [Strike Price], table.cell(colspan: 2,align: center)[$K = 80$], table.cell(colspan: 2,align: center)[$K = 100$], table.cell(colspan: 2,align: center)[$K = 120$],
    [], [Plain MC], [IS], [Plain MC], [IS], [Plain MC], [IS],
    [Theoretical], [0.079477], [], [0.002399], [], [0.000061], [],
    [Estimate], [0.073593], [0.080079], [0.005058], [0.002373], [0.000000], [0.000060],
    [S.E.], [0.008955], [0.000808], [0.002402], [0.000030], [0.000000], [0.000001],
    [R.E.], [12.17%], [1.01%], [47.48%], [1.26%], [NaN%], [1.45%],
    [$hat(mu)$], [], [2.6001], [], [3.6014], [], [4.4569]
    
    )
  ] <european-call-option-results>

  As the strike price $K$ increases, the option becomes increasingly out-of-the-money, making accurate pricing via standard Monte Carlo increasingly difficult due to the rarity of exercise events.
  The importance sampling estimator significantly outperforms the plain Monte Carlo estimator in these scenarios, achieving relative errors of approximately 1% compared to 12.17% and 47.48% for the plain Monte Carlo at $K = 80$ and $K = 100$ respectively.
  For $K = 120$, the plain Monte Carlo estimator fails to produce a meaningful estimate due to the extreme rarity of exercise events, while the importance sampling estimator maintains a low relative error of 1.45%.
  This demonstrates the effectiveness of importance sampling in focusing computational effort on the critical tail regions of the underlying asset's distribution.

  #figure(
    caption: [European Call Option: payoff distribution under original and importance sampling distributions, with $K = 120$ and mode matching.],
  )[
    #image(
      "assets/importance-sampling-european-call-option.png",
      width: 20cm,
    )
  ] <importance-sampling-european-call-option>

  The variance reduction can be understood by examining the payoff distributions under the original and importance sampling distributions. When sampling from the alternative distribution $N(mu^*, 1)$, the samples are concentrated in the region where the option is more likely to be exercised, leading to a higher frequency of non-zero payoffs as shown in @importance-sampling-european-call-option. This results in a more accurate estimate of the option price with significantly reduced variance.
  ],
  [*Applying Control Functionals.* Since control functionals is a post-hoc variance reduction technique, it can be easily combined with importance sampling to further reduce variance.
  For the deep out-of-the-money European call option in the previous example, we apply control functionals on top of the importance sampling estimator using the modified RBF kernel with parameters $eta = 0.1, alpha_1 = 0.0001, alpha_2 = 0.5$. The results are summarized in @european-call-option-cf-is-results.

  #figure(
    caption: [Deep Out-of-the-Money European Call Option: relative error of importance sampling estimator with and without control functionals.],
  )[
    #table(
    columns: 7,
    align: center,

    [Strike Price], table.cell(colspan: 2,align: center)[$K = 80$], table.cell(colspan: 2,align: center)[$K = 100$], table.cell(colspan: 2,align: center)[$K = 120$],
    [], [IS], [IS + CF], [IS], [IS + CF], [IS], [IS + CF],
    [Theoretical], [0.079477], [], [0.002399], [], [0.000061], [],
    [Estimate], [0.080079], [0.079603], [0.002373], [0.002408], [0.000060], [0.000061],
    [S.E.], [0.000808], [0.000111], [0.000030], [0.000006], [0.000001], [0.000000],
    [R.E.], [1.01%], [0.1392%], [1.26%], [0.2391%], [1.45%], [0.3519%],
    [$hat(mu)$], table.cell(colspan: 2)[2.6001],table.cell(colspan: 2)[3.6014],  table.cell(colspan: 2)[4.4569]    
    )

  ] <european-call-option-cf-is-results>

  The addition of control functionals on top of the importance sampling estimator yields further variance reduction across all strike prices.
  ]
)

=== Limitations
Like other variance reduction techniques, the effectiveness of importance sampling hinges on the choice of the alternative distribution $g$.
It may be that the parametric family $cal(G)$ is not flexible enough to approximate the ideal zero-variance distribution well, leading to suboptimal variance reduction.
As an example, consider a bimodal target distribution $h(x) f(x)$, a unimodal $g(x; theta)$ may fail to capture both modes adequately. This could result in unreliable (biased) estimates and unreliable standard errors.

In high-dimensional settings, selecting an appropriate $g$ becomes near impossible due to the curse of dimensionality. The volume of the space increases exponentially with dimension, making it difficult for any fixed parametric family to cover the important regions effectively. This can lead to high variance in the importance weights, causing instability in the estimator.

== Applications in Finance <applications-in-finance>
To illustrate the application of importance sampling in finance, we consider the problem of modelling portfolio credit risk, specifically estimating the probability of large portfolio losses due to defaults.

=== Model Setup
Consider a portfolio of $m$ obligors, each with a probability of default $p_i$ over a fixed time horizon and the loss resulting from default of the $i$-th obligor $c_i$. The loss of the portfolio $L$ can be modelled as follows
$ L  = sum_(i = 1)^m c_i bold(1){X_k > x_k} $
where $(X_1, ..., X_m)$ is a multivariate normal vector such that $PP(X_k > x_k) = p_k$. The $x_k$ are latent variables introduced to capture the dependence structure between obligors. Taking each $X_k$ to be standard normal, we set $x_k = Phi^(-1)(1 - p_k)$ where $Phi$ is the CDF of the standard normal distribution so that the marginal default probabilities are satisfied.

Following the notation of #cite(label("portfolio-credit-risk")), the correlation between obligors is modelled via a multi-factor normal copula:
$ X_k = a_(k 1) Z_1 + a_(k 2) Z_2 + ... + a_(k d) Z_d + b_k epsilon_k $
where $Z_j ~ N(0, 1)$ are common factors affecting all obligors, $epsilon_k ~ N(0, 1)$ are idiosyncratic risk factors specific to each obligor, and $a_(k j), b_k$ are factor loadings satisfying $ sum_(j=1)^d a_(k j)^2 + b_k^2 = 1. $
Let $a_k = (a_(k 1), ..., a_(k d))$ be the factor loading vector for obligor $k$.
Then the correlation between $X_k$ and $X_j$, where $k != j$ is given by $a_k a_j^T$ and the conditional default probability of obligor $k$ given the $Z = (Z_1, ..., Z_d)$ is
$ p_k (Z) = P(X_k > x_k | Z) = Phi((a_k Z + Phi^(-1)(p_k))/(b_k)). $
Using this, we can compute the probability of large portfolio losses $P(L >= x)$ for some threshold $x$ via Monte Carlo simulation by drawing samples of the common factors $Z$ and computing the conditional default probabilities.

=== Example
We consider a portfolio of $m = 100$ obligors in a 10-factor model. The parameters have a form similar to #cite(label("portfolio-credit-risk")) 
$ p_k = 0.01 dot.c (1 + sin(16 pi k \/m)), quad k = 1, ..., m; $
$ c_k = (ceil.l 5k \/ m ceil.r)^2, quad k = 1, ..., m. $
For convenience, the factor loadings $a_(k j)$ are taken to be
$ a_(k j) = 1/sqrt(10)(1 - e^(-j)), quad k = 1, ..., m; #h(0.5em) j = 1, ..., 10 $
We will compute $PP(L > x)$ for the range $x in [0, 200]$ and compare the plain Monte Carlo and importance sampling estimates on $N = 10,000$ samples.
#enum(
  [*Plain Monte Carlo.* For each iteration, we draw sample $Z_i ~ N(0, I_d)$ and compute the latent default variables $X_k$ and portfolio loss $L_i$. The estimate of the tail loss probability is given by the expression
  $ hat(p) = 1/N sum_(i=1)^N bold(1){L_i > x}. $
  These estimates and their 95% confidence intervals are reported in  @portfolio-credit-risk-results-plain-is. We observe that for $x > 80$, the estimates are highly unreliable with wide confidence intervals. In order to accurately estimate such rare-event probabilities, we turn to importance sampling.

  #figure(
    caption: [Portfolio Credit Risk: plain Monte Carlo estimates of $PP(L > x)$.],
  )[
    #image(
      "assets/portfolio-credit-risk-plain-is.png",
      width: 14cm,
    )
  ] <portfolio-credit-risk-results-plain-is>
  ],
  [*Importance Sampling.* The importance sampling approach here is the same as the one described in #cite(label("portfolio-credit-risk")). Let $Y_k$ denote the default indicators $bold(1){X_k > x_k}$. Conditional on the common factors $Z$, the indicators $Y_k$ are independent, so the distribution of $L$ given $Z$ can be shifted to increase the likelihood of $L > x$ by modifying the conditional default probabilities $p_k (Z)$.

  This can be achieved by introducing a parameter $theta$ and defining an exponential twist on the conditional default probabilities as
  $ p_(k, theta) (Z) = (p_k (Z) e^(theta c_k)) / (1 + p_k (Z) ( e^(theta c_k) - 1)). $
  To increase the likelihood of large losses, we increase default probabilities of each obligor by replacing $p_k (Z)$ with $p_(k, theta) (Z)$. Under this exponential twist, the likelihood ratio is given by
  $ product_(k = 1)^m ((p_k (Z))/(p_(k, theta) (Z)))^(Y_k) ((1 - p_k (Z))/(1 - p_(k, theta) (Z)))^(1 - Y_k) = exp(-theta (Z) L + psi(theta, Z)) $
  where $psi(theta, Z) = sum_(k=1)^m log(1 + p_k (Z) (e^(theta c_k) - 1))$. Then the importance sample for $PP(L > x)$ conditional on $Z$ is
  $ bold(1){L_i > x} exp(-theta (Z) L_i + psi(theta, Z)), $
  where the $L_i$ are computed with $p_(k, theta) (Z)$ from $Y_k ~ "Bernoulli"(p_(k, theta) (Z))$.

  The choice of $theta$ depends on $Z$, and is selected to minimize the variance of the importance sampling estimator. This is hard to do directly, but #cite(label("portfolio-credit-risk")) shows that a near-optimal $theta^*(Z)$ can be obtained by solving the equation
  $ partial / (partial theta) psi(theta, Z) = x, $
  in the case where $x >= E[L|Z]$ and setting $theta^*(Z) = 0$ otherwise. 
  
  This method corresponds to minimizing an upper bound on the second moment of the importance sampling estimator #cite(label("portfolio-credit-risk")). 
  
  The results of the importance sampling estimates along with their 95% confidence intervals are reported in @portfolio-credit-risk-results-plain-is.
  We observe that the importance sampling estimates remain reliable even for large thresholds $x approx 200$, demonstrating the effectiveness of importance sampling in estimating rare-event probabilities in portfolio credit risk. 
  ],
  [*Two-step Importance Sampling.* It is worth noting that in situations where the correlation between obligors is high, and the number of obligors $m$ is large, the exponential twist may not be sufficient by itself. To address this, one can employ a two-step importance sampling procedure by also shifting the distribution of the common factors $Z$. This involves introducing another parameter $mu$ and sampling \ $Z ~ N(mu, I_d)$ instead of $N(0, I_d)$. This parameter $mu$ can be selected via an iterative cross-entropy scheme as described in #cite(label("textbook")). 
  This two-step approach can further enhance variance reduction in challenging scenarios with high correlation structures among obligors.
  It can also be noted that when the number of factors $d$ is large, selecting an appropriate shift $mu$ becomes increasingly difficult due to the curse of dimensionality. A full treatment of this two-step importance sampling method is beyond the scope of this report, but interested readers are referred to #cite(label("portfolio-credit-risk")) and #cite(label("textbook")) for more details. 
  ]
)

= Conclusion
This report has explored several variance reduction techniques and their applications in financial engineering. We begain with antithetic sampling and control variates, which are classical methods that leverage negative correlation and known expectations of auxiliary variables to reduce variance. 
Building on these foundations, we delved into control functionals, a more recent development in variance reduction techniques with nice theoretical properties. Control functionals leverage the structure of reproducing kernel Hilbert spaces to construct non-parametric control variates, achieving super-root-$N$ convergence.
The developement of other non-parametric and neural network-based control variates based on similar principles is an active area of research.
Finally, we discussed importance sampling, which, focuses sampling efforts on critical regions of the domain, significantly improving estimation accuracy for rare-event probabilities.
Each technique has its strengths and limitations, and their effectiveness often depends on the specific characteristics of the problem at hand. In practice, these methods can be combined to further enhance variance reduction, as demonstrated with control functionals applied on top of importance sampling.


#pagebreak()
#bibliography(
  "references.bib",
  full: true
)

#pagebreak()
#include("appendix.typ")