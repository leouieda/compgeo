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
    p[1][1] = -sqrt_x
    for n in range(2, n_max + 1):
        p[n] = {}
        for m in range(0, n - 1):
            p[n][m] = (2 * n - 1) / (n - m) * x * p[n - 1][m] - (n + m - 1) / (
                n - m
            ) * p[n - 2][m]
        # Do the diagonal first because we need it for m = n - 1
        p[n][n] = -(2 * n - 1) * sqrt_x * p[n - 1][n - 1]
        m = n - 1
        p[n][m] = (
            -sqrt_x * p[n][m + 1] - sqrt_x * (n + m) * (n - m + 1) * p[n][m - 1]
        ) / (2 * m * x)
    return p


def associated_legendre_functions_derivative(p):
    """
    Calculate d/dx of Plm
    """
    x = p[1][0]
    sqrt_x = np.sqrt(1 - x**2)
    dp = {}
    dp[0] = {}
    dp[0][0] = 0
    dp[1] = {}
    dp[1][0] = 1
    dp[1][1] = x / sqrt_x
    for n in range(2, len(p)):
        dp[n] = {}
        for m in range(0, n):
            dp[n][m] = (n * x * p[n][m] - (n + m) * p[n - 1][m]) / (x**2 - 1)
        m = n
        dp[n][m] = (
            -(n + m) * (n - m + 1) * sqrt_x * p[n][m - 1] - m * x * p[n][m]
        ) / (x**2 - 1)
    return dp


def schmidt_normalization(n_max):
    """
    Calculate the Schmidt normalization factor for geomagnetic field models
    """
    s = {}
    for n in range(n_max + 1):
        s[n] = {}
        for m in range(n + 1):
            delta = 0
            if m == 0:
                delta = 1
            s[n][m] = np.sqrt(
                (2 - delta) * math.factorial(n - m) / math.factorial(n + m)
            )
    return s


def test_legendre_functions():
    "Check if the first few degrees match analytical expressions"
    for angle in np.linspace(0, np.pi, 50):
        x = np.cos(angle)
        p = associated_legendre_functions(x, 4)

        np.testing.assert_allclose(1, p[0][0])

        np.testing.assert_allclose(x, p[1][0])
        np.testing.assert_allclose(-np.sin(angle), p[1][1], atol=1e-10)

        np.testing.assert_allclose(1 / 2 * (3 * x**2 - 1), p[2][0], atol=1e-10)
        np.testing.assert_allclose(-3 * x * np.sin(angle), p[2][1], atol=1e-10)
        np.testing.assert_allclose(3 * (1 - x**2), p[2][2], atol=1e-10)

        np.testing.assert_allclose(1 / 2 * (5 * x**3 - 3 * x), p[3][0], atol=1e-10)
        np.testing.assert_allclose(
            3 / 2 * (1 - 5 * x**2) * np.sin(angle), p[3][1], atol=1e-10
        )
        np.testing.assert_allclose(15 * x * (1 - x**2), p[3][2], atol=1e-10)
        np.testing.assert_allclose(-15 * np.sin(angle) ** 3, p[3][3], atol=1e-10)

        np.testing.assert_allclose(
            1 / 8 * (35 * x**4 - 30 * x**2 + 3), p[4][0], atol=1e-10
        )
        np.testing.assert_allclose(
            -5 / 2 * (7 * x**3 - 3 * x) * np.sin(angle), p[4][1], atol=1e-10
        )
        np.testing.assert_allclose(
            15 / 2 * (7 * x**2 - 1) * (1 - x**2), p[4][2], atol=1e-10
        )
        np.testing.assert_allclose(-105 * x * np.sin(angle) ** 3, p[4][3], atol=1e-10)
        np.testing.assert_allclose(105 * np.sin(angle) ** 4, p[4][4], atol=1e-10)


def test_legendre_functions_derivatives():
    "Check if the first few degrees match analytical expressions"
    for angle in np.linspace(0.1, np.pi - 0.1, 50):
        x = np.cos(angle)
        p = associated_legendre_functions(x, 4)
        dp = associated_legendre_functions_derivative(p)

        np.testing.assert_allclose(0, dp[0][0])

        np.testing.assert_allclose(1, dp[1][0])
        np.testing.assert_allclose(x / np.sqrt(1 - x**2), dp[1][1], atol=1e-10)

        np.testing.assert_allclose(3 * x, dp[2][0], atol=1e-10)
        np.testing.assert_allclose(-3 * np.sqrt(1 - x**2) + 3 * x**2 / np.sqrt(1 - x**2), dp[2][1], atol=1e-10)
        np.testing.assert_allclose(-6 * x, dp[2][2], atol=1e-10)

        np.testing.assert_allclose(15 / 2 * x**2 - 3 / 2, dp[3][0], atol=1e-10)
        np.testing.assert_allclose(
            3 * x * (15 * x**2 - 11) / (2 * np.sqrt(1 - x**2)), dp[3][1], atol=1e-10
        )
        np.testing.assert_allclose(15 - 45 * x**2, dp[3][2], atol=1e-10)
        np.testing.assert_allclose(45 * x * np.sqrt(1 - x**2), dp[3][3], atol=1e-10)


def test_schmidt():
    "Check if the first few degrees match analytical expressions"
    s = schmidt_normalization(3)
    np.testing.assert_allclose(1, s[0][0])
    np.testing.assert_allclose(1, s[1][0])
    np.testing.assert_allclose(1, s[1][1])
    np.testing.assert_allclose(1, s[2][0])
    np.testing.assert_allclose(np.sqrt(1/3), s[2][1])
    np.testing.assert_allclose(np.sqrt(1/12), s[2][2])
