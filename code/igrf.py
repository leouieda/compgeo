"""
Week 6: Calculate the IRGF field at a given point
"""
import datetime
import pytest
import numpy as np
import boule as bl
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

    p = legendre.associated_legendre_functions(np.cos(colatitude), n_max=len(g))
    dp = legendre.associated_legendre_functions_derivative(p)
    s = legendre.schmidt_normalization(n_max=len(g))

    be, bn, br = 0, 0, 0
    for n in g:
        r_frac = (EARTH_RADIUS / radius)**(n + 1)
        r_frac2 = (n + 1) * (EARTH_RADIUS / radius)**(n + 2)
        for m in g[n]:
            cos = np.cos(m * lon_rad)
            sin = np.sin(m * lon_rad)
            if m == 0:
                h[n][m] = 0
            be += r_frac * (-m * g[n][m] * sin + m * h[n][m] * cos) * s[n][m] * p[n][m]
            bn += r_frac * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * dp[n][m]
            br += r_frac2 * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * p[n][m]
    be *= EARTH_RADIUS / (radius * np.sin(colatitude))
    bn *= -EARTH_RADIUS * np.sin(colatitude) / radius

    return be, bn, br


if __name__ == "__main__":
    lon, lat, r = bl.WGS84.geodetic_to_spherical(0, 0, 0)
    be, bn, br = igrf(lon, lat, r, date=datetime.datetime(2020, 1, 1))
    print(bn, be, br, np.sqrt(be**2 + bn**2 + br**2))
    br, bn, be = ppigrf.igrf_gc(r*1e-3, 90 - lat, lon, date=datetime.datetime(2020, 1, 1))
    print(bn, be, br, np.sqrt(be**2 + bn**2 + br**2))
