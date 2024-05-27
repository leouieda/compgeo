"""
Week 8: Make the IGRF code faster
"""
import time
import datetime
import pytest
import numpy as np
import boule as bl
import xarray as xr
import verde as vd
import pygmt

import read_coef
import coef_date
import legendre
import line_profiler



EARTH_RADIUS = 6371.2e3  # m


@line_profiler.profile
def igrf_grid(region, spacing, height, date):
    """
    Make a grid of Bx, By, Bz at a uniform height.
    """
    g_all, h_all, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    g, h = coef_date.coef_date(date, g_all, h_all, g_sv, h_sv, years)

    longitude, latitude = vd.grid_coordinates(region, spacing=spacing, meshgrid=False)

    shape = (latitude.size, longitude.size)
    be = np.zeros(shape)
    bn = np.zeros(shape)
    bu = np.zeros(shape)

    for i, lat in enumerate(latitude):
        for j, lon in enumerate(longitude):
            be_point, bn_point, bu_point = _igrf_internal(lon, lat, height, g, h)
            be[i, j] = be_point
            bn[i, j] = bn_point
            bu[i, j] = bu_point

    dims = ("latitude", "longitude")
    grid = xr.Dataset(
        {
            "be": (dims, be),
            "bn": (dims, bn),
            "bu": (dims, bu),
        },
        coords={"longitude": longitude, "latitude": latitude},
    )
    return grid


@line_profiler.profile
def _igrf_internal(longitude, latitude, height, g, h):
    """
    Internal calculation common to both functions
    """
    longitude, latitude_gc, radius = bl.WGS84.geodetic_to_spherical(longitude, latitude, height)

    colatitude = np.radians(90 - latitude_gc)
    longitude_rad = np.radians(longitude)

    n_max = len(g)
    p = legendre.associated_legendre_functions(np.cos(colatitude), n_max)
    dp = legendre.associated_legendre_functions_derivative(p)
    s = legendre.schmidt_normalization(n_max)

    be, bn_gc, br = 0, 0, 0
    for n in range(1, n_max + 1):
        r_frac = (EARTH_RADIUS / radius)**(n + 2)
        for m in range(0, n + 1):
            cos = np.cos(m * longitude_rad)
            sin = np.sin(m * longitude_rad)
            if m == 0:
                h[n][m] = 0
            be += r_frac * (-m * g[n][m] * sin + m * h[n][m] * cos) * s[n][m] * p[n][m]
            bn_gc += r_frac * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * dp[n][m]
            br += (n + 1) * r_frac * (g[n][m] * cos + h[n][m] * sin) * s[n][m] * p[n][m]
    be *= -1 / np.sin(colatitude)

    # Rotate the vector from geocentric to geodetic
    cos = np.cos(-np.radians(latitude - latitude_gc))
    sin = np.sin(-np.radians(latitude - latitude_gc))
    bn = cos * bn_gc + sin * br
    bu = -sin * bn_gc + cos * br

    return be, bn, bu


@line_profiler.profile
def igrf(longitude, latitude, height, date):
    """
    """
    g_all, h_all, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    g, h = coef_date.coef_date(date, g_all, h_all, g_sv, h_sv, years)

    be, bn, bu = _igrf_internal(longitude, latitude, height, g, h)
    return be, bn, bu


def test_igrf():
    "Check calculated results against those produced by the NOAA website"
    noaa = {
        2018: (22205.3, 3028.0, 45599.4, 50808.9),
        2019: (22193.7, 3060.5, 45678.5, 50876.8),
        2020: (22182.0, 3093.0, 45757.6, 50944.8),
        2021: (22175.3, 3123.3, 45846.0, 51023.1),
        2022: (22168.6, 3153.5, 45934.3, 51101.4),
        2023: (22161.8, 3183.8, 46022.7, 51179.8),
        2024: (22155.1, 3214.1, 46111.0, 51258.2),
    }
    for year in range(2018, 2025):
        date = datetime.datetime(year, 1, 1)
        be, bn, bu = igrf(longitude=45, latitude=45, height=0, date=date)
        amp = np.sqrt(be**2 + bn**2 + bu**2)
        np.testing.assert_allclose(be, noaa[year][1], atol=0.2)
        np.testing.assert_allclose(bn, noaa[year][0], atol=0.2)
        np.testing.assert_allclose(bu, -noaa[year][2], atol=0.2)
        np.testing.assert_allclose(amp, noaa[year][3], atol=0.2)


if __name__ == "__main__":
    start = time.time()
    grid = igrf_grid((0, 360, -90, 90), spacing=3, height=1000, date=datetime.datetime.today())
    print(f"{time.time() - start} s")

    grid["amplitude"] = np.sqrt(grid.be**2 + grid.bn**2 + grid.bu**2)
    fig = pygmt.Figure()
    fig.grdimage(grid.amplitude, projection="W0/20c", cmap="viridis", frame=True)
    fig.colorbar()
    fig.coast(shorelines=True)
    fig.show()
