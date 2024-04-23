"""
Week 2: Creating functions and testing them

Make a function out of the coefficient reading code from last week.
"""
import pathlib
import pytest


def read_gauss_coeffs(path):
    """
    Read Gauss coefficients from the NOAA IGRF data file
    """
    if not pathlib.Path(path).exists():
        raise IOError(f"Gauss coefficient file '{path}' not found.")
    with open(path) as coef_file:
        for i in range(3):
            coef_file.readline()
        line = coef_file.readline()
        years = []
        for year in line.split()[3:-1]:
            years.append(float(year))
        h = {}
        g = {}
        h_sv = {}
        g_sv = {}
        for line in coef_file:
            parts = line.split()
            degree = int(parts[1])
            order = int(parts[2])
            if parts[0] == "g":
                if degree not in g:
                    g[degree] = {}
                    g_sv[degree] = {}
                g[degree][order] = []
                for coef in parts[3:-1]:
                    g[degree][order].append(float(coef))
                g_sv[degree][order] = float(parts[-1])
            if parts[0] == "h":
                if degree not in h:
                    h[degree] = {}
                    h_sv[degree] = {}
                h[degree][order] = []
                for coef in parts[3:-1]:
                    h[degree][order].append(float(coef))
                h_sv[degree][order] = float(parts[-1])
    return g, h, g_sv, h_sv, years


def test_read_all_degrees():
    "Check if all degrees were read"
    g, h, g_sv, h_sv, years = read_gauss_coeffs("igrf13coeffs.txt")
    assert sorted(g.keys()) == list(range(1, 14))
    assert sorted(h.keys()) == list(range(1, 14))
    assert sorted(g_sv.keys()) == list(range(1, 14))
    assert sorted(h_sv.keys()) == list(range(1, 14))


def test_read_all_orders():
    "Check if all orders were read"
    g, h, g_sv, h_sv, years = read_gauss_coeffs("igrf13coeffs.txt")
    for degree in g:
        assert sorted(g[degree].keys()) == list(range(0, degree + 1))
        assert sorted(g_sv[degree].keys()) == list(range(0, degree + 1))
    for degree in h:
        assert sorted(h[degree].keys()) == list(range(1, degree + 1))
        assert sorted(h_sv[degree].keys()) == list(range(1, degree + 1))


def test_years():
    "Check that all years were read"
    g, h, g_sv, h_sv, years = read_gauss_coeffs("igrf13coeffs.txt")
    assert years == list(range(1900, 2025, 5))


def test_file_not_found():
    "Check if it fails when given a bad file name"
    with pytest.raises(IOError):
        read_gauss_coeffs("bla.slkdjsldjh")
