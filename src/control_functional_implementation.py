import numpy as np
from scipy import linalg


def modified_rbf_kernel(X, Y, length_scale=1.0, decay_scale=1.0):
    """Compute the modified RBF kernel matrix between two sets of points."""
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    X_norm2 = np.sum(X**2, axis=1).reshape(-1, 1)
    Y_norm2 = np.sum(Y**2, axis=1).reshape(1, -1)
    sqdist = X_norm2 + Y_norm2 - 2.0 * (X @ Y.T)
    sqdist = np.maximum(sqdist, 0.0)
    rbf = np.exp(-0.5 * sqdist / (length_scale**2))
    decay = 1.0 / (1.0 + decay_scale * X_norm2 + decay_scale * Y_norm2)
    return rbf * decay


def modified_stein_rbf_kernel(X, Y, length_scale=1.0, decay_scale=1e-3):
    """Compute the Stein kernel matrix using the modified RBF kernel."""
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)

    n, p = X.shape
    m = Y.shape[0]

    product = X @ Y.T
    X_norm2 = np.sum(X**2, axis=1).reshape(-1, 1)
    Y_norm2 = np.sum(Y**2, axis=1).reshape(1, -1)
    sqdist = X_norm2 + Y_norm2 - 2.0 * (X @ Y.T)

    s2 = length_scale**2

    rbf = np.exp(-0.5 * sqdist / s2)
    decay = 1.0 / (1.0 + decay_scale * X_norm2 + decay_scale * Y_norm2)

    factor = (
        p / s2
        - sqdist / (s2 * s2)
        - sqdist / s2 * (1.0 + 2.0 * decay_scale * decay)
        + 4.0 * decay_scale * decay * product
        + product
    )

    return rbf * decay * factor


class ControlFunctional:
    """A Control Functional model for variance reduction in Monte Carlo integration."""
    def __init__(self, X, Y, kernel_fn=modified_stein_rbf_kernel, eta=0.10, length_scale=1.0, decay_scale=1e-3):
        self.kernel_fn = kernel_fn
        self.length_scale = length_scale
        self.decay_scale = decay_scale
        
        # Split into training and test sets
        self.n = X.shape[0]
        self.m = int(self.n * eta)

        perm = np.random.permutation(self.n)
        X = X[perm]
        Y = Y[perm]

        self.X_train = X[:self.m]
        self.Y_train = Y[:self.m]
        self.X_test = X[self.m:]
        self.Y_test = Y[self.m:]

        # Placeholders
        self.coeffs = None
        self.const_fn_proj = None

        # Estimators
        self.unbiased_estimator = None
        self.simplified_estimator = None
        self.std_err = None
        self.re_err = None

        # Predictions and Residuals
        self.pred = None
        self.residuals = None

    def fit(self, alpha=1e-3):
        """Fit a function to the training data using kernel ridge regression."""
        self.alpha = alpha
        lam = alpha * np.sqrt(self.m)  # regularization parameter

        K = self.kernel_fn(self.X_train, self.X_train, length_scale=self.length_scale, decay_scale=self.decay_scale)
        Klam = K + lam * np.eye(self.m)

        # Cholesky decomposition and solve for coefficients
        c, lower = linalg.cho_factor(Klam, lower=True)
        self.coeffs = linalg.cho_solve((c, lower), self.Y_train)

        # Cholesky decomposition and solve for normalization constant
        ones = np.ones((self.m, 1))
        self.const_fn_proj = linalg.cho_solve((c, lower), ones)

    def predict(self):
        """Predict on the test set and compute estimators."""
        self.simplified_estimator = np.sum(self.coeffs) / (1.0 + np.sum(self.const_fn_proj))
        if self.n == self.m:
            return

        K_test = self.kernel_fn(self.X_test, self.X_train, length_scale=self.length_scale, decay_scale=self.decay_scale)

        ones = np.ones((self.n - self.m, 1))
        self.pred = K_test @ self.coeffs + (ones - K_test @ self.const_fn_proj) * self.simplified_estimator

        self.unbiased_estimator = np.mean(self.Y_test - self.pred) + self.simplified_estimator
        self.residuals = self.Y_test - self.pred

        # Standard error and relative error
        self.std_err = np.std(self.residuals, ddof=1) / np.sqrt(self.n - self.m)
        eps = 1e-12
        self.re_err = self.std_err / (np.abs(self.unbiased_estimator) + eps)

        return self.pred

    def summary(self):
        """Print a summary of the results."""
        if self.n != self.m:
            print(f"Unbiased Est.: {self.unbiased_estimator:.6f}")
            print(f"Simplified Est.: {self.simplified_estimator:.6f}")
            print(f"SE: {self.std_err:.6f}")
            print(f"RE: {self.re_err:.4%}")
        else:
            print(f"Simplified Est.: {self.simplified_estimator:.6f}")

    def diagnostics(self):
        # Check matrix condition number
        K = self.kernel_fn(self.X_train, self.X_train, length_scale=self.length_scale, decay_scale=self.decay_scale)
        cond_number = np.linalg.cond(K)
        print(f"Condition Number of Kernel Matrix: {cond_number:.2e}")

        # Check K + lambda I condition number
        lam = self.alpha * np.sqrt(self.m)
        Klam = K + lam * np.eye(self.m)
        cond_number_lam = np.linalg.cond(Klam)
        print(f"Condition Number of Regularized Kernel Matrix: {cond_number_lam:.2e}")

        # Coefficients magnitude
        coeffs_norm = np.linalg.norm(self.coeffs)
        print(f"Norm of Coefficients: {coeffs_norm:.2e}")