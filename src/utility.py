import numpy as np

from scipy.stats import multivariate_normal


def bisection_method(h, x1, x2, eps=1e-5):
    """Finds a root of a function h within a given interval using the bisection method.

    Assumes that the function h is continuous and that
    h(x1) and h(x2) have opposite signs (i.e., the root is bracketed
    in the interval [x1, x2]).
    """
    # Check that h(x1) and h(x2) have opposite signs.
    assert h(x1) * h(x2) < 0

    while abs(x1 - x2) >= eps:
        x = (x1 + x2) / 2
        if x == 0:
            x1 = x
            x2 = x
        else:
            if h(x) * h(x1) > 0:
                x1 = x
            else:
                x2 = x
    
    return x1

def secant_method(h, x1, x2, eps=1e-5, max_iter=1000):
    """Finds a root of a function h using the secant method.
    """
    for _ in range(max_iter):
        if abs(h(x2)) < eps: return x2
        if h(x2) == h(x1): break
        x3 = x2 - h(x2) * (x2 - x1) / (h(x2) - h(x1))
        x1, x2 = x2, x3
    return x2


def norm_cross_entropy(h, m, N):
    """Estimate the minimizer of the KL-divergence between a normal distribution
    and a reweighted distribution defined by function h.

    Given a measurable function h: R^m -> R, this function generates N samples
    from the m-dimensional standard normal distribution, reweights them using h,
    and computes the empirical minimizer of the KL-divergence (from Lemma 7.1 in the
    referenced context).
    """
    X = multivariate_normal.rvs(mean=np.zeros(m), cov=np.eye(m), size=N)
    if m == 1:
        X = X.reshape(-1, 1)
    hX = np.apply_along_axis(h, 1, X)
    theta_hat = np.sum(hX.reshape(-1, 1) * X, axis=0) / np.sum(hX)
    if m ==1 :
        return theta_hat[0]
    return theta_hat

def iterative_norm_cross_entropy(h, m, N, theta):
    Y = multivariate_normal.rvs(mean=theta, cov=np.eye(m), size=N)
    if m == 1:
        Y = Y.reshape(-1, 1)
    hY = np.apply_along_axis(h, 1, Y)
    e_term = np.exp(-np.dot(Y, theta))
    hYe = hY * e_term
    theta_hat = np.sum(hYe.reshape(-1, 1) * Y  , axis=0) / np.sum(hYe)
    if m == 1 :
        return theta_hat[0]
    return theta_hat

def control_variate_coeffient(X, Y):
    """Compute the optimal control variate coefficient b* given two arrays X and Y.
    """
    cov_matrix = np.cov(X, Y)
    b_star = cov_matrix[0, 1] / cov_matrix[1, 1]
    return b_star