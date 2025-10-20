import jax
import jax.numpy as jnp

@jax.jit
def rbf_kernel(X, Y, length_scale=1.0):
    """Compute the RBF kernel matrix between two sets of points."""
    X = jnp.atleast_2d(X)
    Y = jnp.atleast_2d(Y)
    X_norm2 = jnp.sum(X**2, axis=1).reshape(-1, 1)
    Y_norm2 = jnp.sum(Y**2, axis=1).reshape(1, -1)
    sqdist = X_norm2 + Y_norm2 - 2.0 * (X @ Y.T)
    sqdist = jnp.maximum(sqdist, 0.0)
    return jnp.exp(-0.5 * sqdist / (length_scale**2))

@jax.jit
def rbf_stein_kernel(X, Y, length_scale=1.0):
    """Compute the Stein kernel matrix using the RBF kernel."""
    X = jnp.atleast_2d(X)
    Y = jnp.atleast_2d(Y)
    product = X @ Y.T
    X_norm2 = jnp.sum(X**2, axis=1).reshape(-1, 1)
    Y_norm2 = jnp.sum(Y**2, axis=1).reshape(1, -1)
    sqdist = X_norm2 + Y_norm2 - 2.0 * product
    sqdist = jnp.maximum(sqdist, 0.0)
    rbf = jnp.exp(-0.5 * sqdist / (length_scale**2))
    factor = (1.0 / length_scale**2) - sqdist * (1.0 / length_scale**4 + 1.0 / length_scale**2) + product
    return rbf * factor

class ControlFunctional:
    """A Control Functional model for variance reduction in Monte Carlo integration."""
    def __init__(self, X, Y, key, kernel_fn=rbf_stein_kernel, r=0.10, length_scale=1.0):
        self.kernel_fn = kernel_fn
        self.length_scale = length_scale
        
        # Split into training and test sets
        self.n = X.shape[0]
        self.m = int(self.n * r)

        perm = jax.random.permutation(key, self.n)
        X = X[perm]
        Y = Y[perm]

        self.X_train = X[:self.m]
        self.Y_train = Y[:self.m]
        self.X_test = X[self.m:]
        self.Y_test = Y[self.m:]

        # Placeholder for fitted coefficients
        self.coeffs = None
        self.const_fn_proj = None

        # Estimators
        self.unbiased_estimator = None
        self.simplified_estimator = None
        
        self.std_err = None
        self.re_err = None

        # Predictions and Residuals
        self.pred = None
        self.resid = None
    
    def fit(self, alpha=1e-3):
        """Fit a function to the training data using kernel ridge regression."""
        lam = alpha * jnp.sqrt(self.m) # regularization parameter

        K = self.kernel_fn(self.X_train, self.X_train, length_scale=self.length_scale)
        Klam = K + lam * jnp.eye(self.m)

        # Cholesky decomposition and solve for coefficients
        L = jnp.linalg.cholesky(Klam)
        z = jax.scipy.linalg.solve_triangular(L, self.Y_train.reshape(-1, 1), lower=True)
        self.coeffs = jax.scipy.linalg.solve_triangular(L.T, z, lower=False)

        # Cholesky decomposition and solve for normalization constant
        ones = jnp.ones((self.m, 1))
        z_const = jax.scipy.linalg.solve_triangular(L, ones, lower=True)
        self.const_fn_proj = jax.scipy.linalg.solve_triangular(L.T, z_const, lower=False)


    def predict(self):
        """Predict on the test set and compute estimators."""
        self.simplified_estimator = jnp.sum(self.coeffs) / (1.0 + jnp.sum(self.const_fn_proj))
        if self.n == self.m:
            return

        K_test = self.kernel_fn(self.X_test, self.X_train, length_scale=self.length_scale)

        ones = jnp.ones((self.n - self.m, 1))
        self.pred = K_test @ self.coeffs + (ones - K_test @ self.const_fn_proj) * self.simplified_estimator

        self.unbiased_estimator = jnp.mean(self.Y_test - self.pred) + self.simplified_estimator
        self.residuals = self.Y_test - self.pred

        # Standard error and relative error
        self.std_err = jnp.sqrt(jnp.mean(self.residuals**2) / (self.n - self.m))
        self.re_err = self.std_err / jnp.abs(self.unbiased_estimator)


    def summary(self):
        """Print a summary of the results."""
        if self.n != self.m:
            print(f"Unbiased Est.: {self.unbiased_estimator:.6f}")
            print(f"Simplified Est.: {self.simplified_estimator:.6f}")
            print(f"SE: {self.std_err:.6f}")
            print(f"RE: {self.re_err:.4%}")
        else:
            print(f"Simplified Est.: {self.simplified_estimator:.6f}")
    