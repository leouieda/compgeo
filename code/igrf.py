"""
Week 6: Calculate the IRGF field at a given point
"""
import datetime
import pytest
import numpy as np
import boule as bl
import scipy
import ppigrf
import read_coef
import coef_date
import legendre


EARTH_RADIUS = 6371.2e3  # m


def igrf(longitude, latitude, radius, date):
    """
    """
    g_all, h_all, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    g, h = coef_date.coef_date(date, g_all, h_all, g_sv, h_sv, years)

    colatitude = np.radians(90 - latitude)
    lon_rad = np.radians(longitude)

    n_max = len(g)
    p = legendre.associated_legendre_functions(np.cos(colatitude), n_max)
    dp = legendre.associated_legendre_functions_derivative(p)
    s = legendre.schmidt_normalization(n_max)

    be, bn, br = 0, 0, 0
    for n in range(1, n_max + 1):
        r_frac = (EARTH_RADIUS / radius)**(n + 2)
        for m in range(0, n + 1):
            cos = np.cos(m * lon_rad)
            sin = np.sin(m * lon_rad)
            if m == 0:
                h[n][m] = 0
            be += r_frac * (-m * g[n][m] * sin + m * h[n][m] * cos) * s[n][m] * p[n][m]
            bn += r_frac * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * dp[n][m]
            br += (n + 1) * r_frac * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * p[n][m]
    be *= -1 / np.sin(colatitude)
    return be, bn, br


if __name__ == "__main__":
    lon, lat, r = bl.WGS84.geodetic_to_spherical(45, 45, 0)
    date = datetime.datetime(2000, 1, 1)
    be, bn, br = igrf(lon, lat, r, date)
    print(bn, be, br, np.sqrt(be**2 + bn**2 + br**2))
    br, bn, be = ppigrf.igrf_gc(r*1e-3, 90 - lat, lon, date)
    print(bn[0, 0], be[0, 0], br[0, 0], np.sqrt(be**2 + bn**2 + br**2)[0, 0])
