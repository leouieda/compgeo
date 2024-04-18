"""
Function to read the Gauss coefficients along with some tests.
"""
import pathlib
import pytest


def read_gauss_coeffs(path):
    """
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
    # Add the zero order h coefficients
    for degree in h:
        h[degree][0] = [0] * len(years)
        h_sv[degree][0] = 0
    return g, h, g_sv, h_sv, years


def test_read_all_degrees():
    "Check if all degrees were read"
    g, h, g_sv, h_sv, years = read_gauss_coeffs("igrf13coeffs.txt")
    assert sorted(g.keys()) == list(range(1, 14))
    assert sorted(h.keys()) == list(range(1, 14))
    assert sorted(g_sv.keys()) == list(range(1, 14))
    assert sorted(h_sv.keys()) == list(range(1, 14))


def test_file_not_found():
    "Check if it fails when given a bad file name"
    with pytest.raises(IOError):
        read_gauss_coeffs("bla.slkdjsldjh")
