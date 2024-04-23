"""
Week 3: Calculate the dipole moment and the pole locations
"""
import pytest
import pygmt
import numpy as np
import read_coef


VACUUM_PERMEABILITY = 1.2566370621219e-6  # N/A²
EARTH_RADIUS = 6.3781e6  # m


def coef_year(year, g, h, g_sv, h_sv, years):
    """
    Get the coefficients for a specific year
    """
    if year not in years:
        raise ValueError(f"Invalid year '{year}'. Must be one of {years}")
    index = years.index(year)
    g_year = {}
    h_year = {}
    for n in g:
        g_year[n] = {}
        for m in g[n]:
            g_year[n][m] = g[n][m][index]
    for n in h:
        h_year[n] = {}
        for m in h[n]:
            h_year[n][m] = h[n][m][index]
    return g_year, h_year


def dipole_moment(g, h):
    """
    Calculate the dipole moment given the Gauss coefficients for a given date
    """
    mx = 4 * np.pi / VACUUM_PERMEABILITY * EARTH_RADIUS**3 * g[1][1]
    my = 4 * np.pi / VACUUM_PERMEABILITY * EARTH_RADIUS**3 * h[1][1]
    mz = 4 * np.pi / VACUUM_PERMEABILITY * EARTH_RADIUS**3 * g[1][0]
    return mx, my, mz


def to_spherical(x, y, z):
    """
    Convert geocentric Cartesian coordinates to spherical
    """
    amplitude = np.sqrt(x**2 + y**2 + z**2)
    latitude = np.degrees(np.arcsin(z / amplitude))
    longitude = np.degrees(np.arctan2(y, x))
    return amplitude, longitude, latitude


def test_invalid_year():
    "Check if raises an error when year is invalid"
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    with pytest.raises(ValueError):
        coef_year(1902, g, h, g_sv, h_sv, years)


def test_g_year():
    "Check if the right values were returned for g"
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    g_year, h_year = coef_year(1905, g, h, g_sv, h_sv, years)
    assert g_year[1][0] == -31464
    assert g_year[7][3] == 33
    g_year, h_year = coef_year(2020, g, h, g_sv, h_sv, years)
    assert g_year[1][0] == -29404.8
    assert g_year[6][3] == -121.5


def test_h_year():
    "Check if the right values were returned for h"
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    g_year, h_year = coef_year(1905, g, h, g_sv, h_sv, years)
    for n in h_year:
        assert 0 not in h_year[n]
    assert h_year[1][1] == 5909
    assert h_year[7][3] == -11
    g_year, h_year = coef_year(2020, g, h, g_sv, h_sv, years)
    for n in h_year:
        assert 0 not in h_year[n]
    assert h_year[1][1] == 4652.5
    assert h_year[6][3] == 52.8


if __name__ == "__main__":
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    lons = []
    lats = []
    for year in years:
        g_year, h_year = coef_year(year, g, h, g_sv, h_sv, years)
        mx, my, mz = dipole_moment(g_year, h_year)
        amp, lon, lat = to_spherical(mx, my, mz)
        lons.append(lon)
        lats.append(lat)

    fig = pygmt.Figure()
    fig.coast(
        region=[0, 360, -90, -60],
        projection="S0/-90/20c",
        land="lightgray",
        frame=True,
    )
    fig.plot(x=lons, y=lats, style="c0.1c", fill="black")
    fig.show()
