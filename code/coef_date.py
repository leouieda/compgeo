"""
Week 4: Interpolating the coefficients for a given date
"""
import datetime
import pytest
import pygmt
import numpy as np
import read_coef
import dipole_moment


def coef_date(date, g, h, g_sv, h_sv, years):
    """
    Get the coefficients for a specific year
    """
    max_date = datetime.datetime(year=int(years[-1]) + 6, month=1, day=1)
    min_date = datetime.datetime(year=int(years[0]), month=1, day=1)
    if date >= max_date or date < min_date:
        raise ValueError(
            f"Invalid date '{str(date)}'. Must be < {max_date} and >= {min_date}",
        )

    index = int((date.year - years[0]) // 5)
    timedelta = date - datetime.datetime(int(years[index]), 1, 1)
    epoch = (
        datetime.datetime(year=int(years[index]) + 5, month=1, day=1)
        - datetime.datetime(year=int(years[index]), month=1, day=1)
    )
    year_in_seconds = (365.25 * 24 * 60 * 60)

    g_year = {}
    h_year = {}
    for n in g:
        g_year[n] = {}
        for m in g[n]:
            if index < len(years) - 1:
                secular_variation = (g[n][m][index + 1] - g[n][m][index]) / epoch.total_seconds()
            else:
                secular_variation = g_sv[n][m] / year_in_seconds
            g_year[n][m] = g[n][m][index] + timedelta.total_seconds() * secular_variation
    for n in h:
        h_year[n] = {}
        for m in h[n]:
            if index < len(years) - 1:
                secular_variation = (h[n][m][index + 1] - h[n][m][index]) / epoch.total_seconds()
            else:
                secular_variation = h_sv[n][m] / year_in_seconds
            h_year[n][m] = h[n][m][index] + timedelta.total_seconds() * secular_variation
    return g_year, h_year


def test_invalid_year():
    "Check if raises an error when year is invalid"
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    with pytest.raises(ValueError):
        coef_date(datetime.datetime.fromisoformat("1899-12-31T23:59:59"), g, h, g_sv, h_sv, years)
    with pytest.raises(ValueError):
        coef_date(datetime.datetime.fromisoformat("2026-01-01T00:00:00"), g, h, g_sv, h_sv, years)


def test_interpolation():
    "Check that calculating on or close to an epoch yields the right coefficients"
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    for i, year in enumerate(years):
        g_date, h_date = coef_date(
            datetime.datetime(int(year), 1, 1, 0, 1, 0), g, h, g_sv, h_sv, years,
        )
        for n in g_date:
            for m in g_date[n]:
                np.testing.assert_allclose(g_date[n][m], g[n][m][i], atol=0.001)
        for n in h_date:
            for m in h_date[n]:
                np.testing.assert_allclose(h_date[n][m], h[n][m][i], atol=0.001)


if __name__ == "__main__":
    g, h, g_sv, h_sv, years = read_coef.read_gauss_coeffs("igrf13coeffs.txt")
    fig = pygmt.Figure()
    fig.coast(
        region=[60, 150, -85, -65],
        projection="S60/-90/20c",
        land="lightgray",
        frame=True,
    )
    for year in range(1900, 2020):
        for month in range(1, 13):
            g_year, h_year = coef_date(datetime.datetime(int(year), month, day=1), g, h, g_sv, h_sv, years)
            mx, my, mz = dipole_moment.dipole_moment(g_year, h_year)
            amp, lon, lat = dipole_moment.to_spherical(mx, my, mz)
            fig.plot(x=lon, y=lat, style="c0.1c", fill="black")
    for year in range(2020, 2025):
        for month in range(1, 13):
            g_year, h_year = coef_date(datetime.datetime(int(year), month, day=1), g, h, g_sv, h_sv, years)
            mx, my, mz = dipole_moment.dipole_moment(g_year, h_year)
            amp, lon, lat = dipole_moment.to_spherical(mx, my, mz)
            fig.plot(x=lon, y=lat, style="c0.1c", fill="blue")
    fig.show()
