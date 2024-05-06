import numpy as np


def associated_legendre_functions(x, n_max):
    """
    """
    sqrt_x = np.sqrt(1 - x**2)

    p = {}
    p[0] = {}
    p[0][0] = 1
    p[1] = {}
    p[1][0] = x
    p[1][1] = -sqrt_x

    for n in range(2, n_max + 1):
        p[n] = {}
        for m in range(0, n + 1):
            if n == m:
                p[n][m] = -(2 * n - 1) * sqrt_x * p[n - 1][m - 1]
            elif m == n - 1:
                p[n][m] = (2 * n - 1) / (n - m) * x * p[n - 1][m]
            else:
                p[n][m] = (
                    (2 * n - 1) / (n - m) * x * p[n - 1][m]
                    - (n + m - 1) / (n - m) * p[n - 2][m]
                )
    return p



