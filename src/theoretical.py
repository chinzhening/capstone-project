import numpy as np

from scipy.stats import norm

def european_call(S_0, K, r, T, sig):
    """
    Price of a European call option under the Black-Scholes model.

    Parameters
    ----------
    S_0 : float
        Current stock price.
    K : float
        Strike price of the option.
    r : float
        Risk-free interest rate (annualized, continuously compounded).
    T : float
        Time to maturity in years.
    sig : float
        Volatility of the underlying asset (annualized standard deviation).

    Returns
    -------
    float
        Present value of the European call option.

    Notes
    -----
    The formula used is the Black-Scholes closed-form solution:
        C = S_0 * N(d1) - K * exp(-r T) * N(d2),
    where
        d1 = [ln(S_0 / K) + (r + 0.5 sig^2) T] / (sig sqrt(T)),
        d2 = d1 - sig sqrt(T).
    """
    d = (np.log(S_0 / K) + (r + 0.5 * sig**2) * T) / (sig * np.sqrt(T))
    return S_0 * norm.cdf(d) - K * np.exp(-r * T) * norm.cdf(d - sig * np.sqrt(T))


# Example 2.3
def binary_call(S_0, K, r, T, sig):
    """
    Price of a European cash-or-nothing binary call option under the Black-Scholes model.

    Parameters
    ----------
    S_0 : float
        Current stock price.
    K : float
        Strike price of the option.
    r : float
        Risk-free interest rate (annualized, continuously compounded).
    T : float
        Time to maturity in years.
    sig : float
        Volatility of the underlying asset (annualized standard deviation).

    Returns
    -------
    float
        Present value of the binary call option, which pays 1 unit of currency
        at maturity if S_T > K, and 0 otherwise.

    Notes
    -----
    The formula used is:
        BinaryCall = exp(-r T) * N(-b),
    where
        b = [ln(K / S_0) - (r - 0.5 sig^2) T] / (sig sqrt(T)).
    """
    b = (np.log(K / S_0) - (r - sig**2 / 2) * T) / (sig * np.sqrt(T))
    return np.exp(-r * T) * norm.cdf(-b)

# Example 2.5
def discretely_monitored_average_price_call_option(S_0, K, r, T, sig, m):
    """
    Price of a discretely monitored average price call option under the Black-Scholes model.

    Parameters
    ----------
    S_0 : float
        Current stock price.
    K : float
        Strike price of the option.
    r : float
        Risk-free interest rate (annualized, continuously compounded).
    T : float
        Time to maturity in years.
    sig : float
        Volatility of the underlying asset (annualized standard deviation).
    m : int
        Number of monitoring dates (including maturity).
    Returns
    -------
    float
        Present value of the average price call option. The payoff at maturity is
        max( (1/m) * sum_{i=1}^m S_{t_i} - K, 0 ), where t_i are the monitoring dates.

    Notes
    -----
    The formula used is:
        C = exp(-r T) * [ exp(mu_bar + 0.5 sig_bar^2) * N(sqrt(sig_bar^2) - theta) - K * N(-theta) ],
    where
        t_bar = T * (m + 1) / (2 * m),
        mu_bar = ln(S_0) + (r - 0.5 sig^2) * t_bar,
        sig_bar^2 = sig^2 * (m + 1)(2m + 1) / (6 m^2) * T,
        theta = [ln(K) - mu_bar] / sqrt(sig_bar^2).
    """
    t_bar = T * (m + 1) / (2 * m)
    mu_bar = np.log(S_0) + (r - 0.5 * sig**2) * t_bar
    sig_bar2 = sig**2 * (m + 1) * (2 * m + 1) / (6 * m**2) * T
    theta = (np.log(K) - mu_bar) / np.sqrt(sig_bar2)
    payoff = np.exp(mu_bar + 0.5 * sig_bar2) * norm.cdf(np.sqrt(sig_bar2) - theta) - K * norm.cdf(-theta)
    return np.exp(-r * T) * payoff

# Example 2.6
def continuous_lookback_call(S_0, K, r, T, sig):
    """
    Price of a continuous lookback call option under the Black-Scholes model.

    Parameters
    ----------
    S_0 : float
        Current stock price.
    K : float
        Strike price of the option.
    r : float
        Risk-free interest rate (annualized, continuously compounded).
    T : float
        Time to maturity in years.
    sig : float
        Volatility of the underlying asset (annualized standard deviation).
    Returns
    -------
    float
        Present value of the continuous lookback call option. The payoff at maturity is
        max_{0 <= t <= T} S_t - K.

    """
    N = norm.cdf
    d = (np.log(S_0 / K) + (r + 0.5 * sig**2) * T) / (sig * np.sqrt(T))
    d1 = (r + 0.5 * sig**2) * np.sqrt(T) / sig
    p = sig**2 * 0.5 / r

    term1 = np.exp(r * T) * N(d) - (S_0 / K)**(-1/p) * N(d - 2*r/(sig * np.sqrt(T)))
    term2 = np.exp(r * T) * N(d1) -  N(d1 - 2*r/(sig * np.sqrt(T)))

    if K >= S_0:
        return european_call(S_0, K, r, T, sig) + np.exp(-r*T) * p * S_0 * term1
    else:
        return S_0 * N(d1) + np.exp(-r*T) * (-S_0 * N(d1 - sig * np.sqrt(T)) + S_0 - K + p * S_0 * term2)

def butterfly_spread(S_0, K, r, T, sig):
    K1, K2, K3 = K
    p1 = european_call(S_0, K1, r, T, sig)
    p2 = european_call(S_0, K2, r, T, sig)
    p3 = european_call(S_0, K3, r, T, sig)
    return p1 + p3 - 2 * p2

def straddle_option(S_0, K, r, T, sig):
    call_price = european_call(S_0, K, r, T, sig)
    put_price = call_price + K * np.exp(-r * T) - S_0
    return call_price + put_price