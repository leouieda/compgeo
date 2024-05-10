"""
Week 5: Calculate associated Legendre functions, their derivatives, and the
Schmidt normalization factor.
"""
import numpy as np
import math


def associated_legendre_functions(x, n_max):
    """
    Calculate Plm for x in range [-1, 1]
    """
    sqrt_x = np.sqrt(1 - x**2)
    p = {}
    p[0] = {}
    p[0][0] = 1
    p[1] = {}
    p[1][0] = x
    p[1][1] = sqrt_x
    for n in range(2, n_max + 1):
        p[n] = {}
        for m in range(0, n - 1):
            p[n][m] = (2 * n - 1) / (n - m) * x * p[n - 1][m] - (n + m - 1) / (
                n - m
            ) * p[n - 2][m]
        m = n - 1
        p[n][m] = (2 * n - 1) / (n - m) * x * p[n - 1][m]
        p[n][n] = (2 * n - 1) * sqrt_x * p[n - 1][n - 1]
    return p


def associated_legendre_functions_derivative(p):
    """
    Calculate d/dtheta of Plm
    """
    x = p[1][0]
    dp = {}
    dp[0] = {}
    dp[0][0] = 0
    for n in range(1, len(p)):
        dp[n] = {}
        dp[n][0] = -p[n][1]
        for m in range(1, n):
            dp[n][m] = 0.5 * ((n + m) * (n - m + 1) * p[n][m - 1] - p[n][m + 1])
        m = n
        dp[n][m] = 0.5 * (n + m) * (n - m + 1) * p[n][m - 1]
    return dp


def schmidt_normalization(n_max):
    """
    Calculate the Schmidt normalization factor for geomagnetic field models
    """
    s = {}
    for n in range(n_max + 1):
        s[n] = {}
        for m in range(n + 1):
            if m == 0:
                s[n][m] = 1
            else:
                s[n][m] = np.sqrt(
                    2 * math.factorial(n - m) / math.factorial(n + m)
                )
    return s


def test_legendre_functions():
    "Check if the first few degrees match analytical expressions"
    for angle in np.linspace(0, np.pi, 50):
        x = np.cos(angle)
        p = associated_legendre_functions(x, 4)

        np.testing.assert_allclose(1, p[0][0])

        np.testing.assert_allclose(x, p[1][0])
        np.testing.assert_allclose(np.sin(angle), p[1][1], atol=1e-10)

        np.testing.assert_allclose(1 / 2 * (3 * x**2 - 1), p[2][0], atol=1e-10)
        np.testing.assert_allclose(3 * x * np.sin(angle), p[2][1], atol=1e-10)
        np.testing.assert_allclose(3 * (1 - x**2), p[2][2], atol=1e-10)

        np.testing.assert_allclose(1 / 2 * (5 * x**3 - 3 * x), p[3][0], atol=1e-10)
        np.testing.assert_allclose(
            -3 / 2 * (1 - 5 * x**2) * np.sin(angle), p[3][1], atol=1e-10
        )
        np.testing.assert_allclose(15 * x * (1 - x**2), p[3][2], atol=1e-10)
        np.testing.assert_allclose(15 * np.sin(angle) ** 3, p[3][3], atol=1e-10)

        np.testing.assert_allclose(
            1 / 8 * (35 * x**4 - 30 * x**2 + 3), p[4][0], atol=1e-10
        )
        np.testing.assert_allclose(
            5 / 2 * (7 * x**3 - 3 * x) * np.sin(angle), p[4][1], atol=1e-10
        )
        np.testing.assert_allclose(
            15 / 2 * (7 * x**2 - 1) * (1 - x**2), p[4][2], atol=1e-10
        )
        np.testing.assert_allclose(105 * x * np.sin(angle) ** 3, p[4][3], atol=1e-10)
        np.testing.assert_allclose(105 * np.sin(angle) ** 4, p[4][4], atol=1e-10)


def test_legendre_functions_derivatives():
    "Check if the first few degrees match analytical expressions"
    for angle in np.linspace(0, np.pi, 50):
        x = np.cos(angle)
        cos = np.cos(angle)
        sin = np.sin(angle)
        p = associated_legendre_functions(x, 4)
        dp = associated_legendre_functions_derivative(p)

        np.testing.assert_allclose(0, dp[0][0])

        np.testing.assert_allclose(-sin, dp[1][0], atol=1e-10)
        np.testing.assert_allclose(cos, dp[1][1], atol=1e-10)

        np.testing.assert_allclose(-3 * cos * sin, dp[2][0], atol=1e-10)
        np.testing.assert_allclose(3 * (cos**2 - sin**2), dp[2][1], atol=1e-10)
        np.testing.assert_allclose(6 * sin * cos, dp[2][2], atol=1e-10)

        np.testing.assert_allclose((3 * sin - 15 * cos**2 * sin) / 2, dp[3][0], atol=1e-10)
        np.testing.assert_allclose(
            3 / 2 * (5 * cos**3 - 10 * cos * sin**2 - cos), dp[3][1], atol=1e-10
        )
        np.testing.assert_allclose(30 * cos**2 * sin - 15 * sin**3, dp[3][2], atol=1e-10)
        np.testing.assert_allclose(45 * sin**2 * cos, dp[3][3], atol=1e-10)


def test_schmidt():
    "Check if the first few degrees match analytical expressions"
    s = schmidt_normalization(3)
    np.testing.assert_allclose(1, s[0][0])
    np.testing.assert_allclose(1, s[1][0])
    np.testing.assert_allclose(1, s[1][1])
    np.testing.assert_allclose(1, s[2][0])
    np.testing.assert_allclose(np.sqrt(1/3), s[2][1])
    np.testing.assert_allclose(np.sqrt(1/12), s[2][2])


if __name__ == "__main__":
    import scipy.special
    n = 13
    x = 0.5
    p, dp = scipy.special.lpmn(n, n, x)
    print(dp.T[-1])
    p = associated_legendre_functions(x, n)
    dp = associated_legendre_functions_derivative(p)
    print(dp[13])
