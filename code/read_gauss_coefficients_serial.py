"""
Week 1: Read the Gauss coefficients from the USGS file. Don't use a function.
"""

h = {}
g = {}
h_sv = {}
g_sv = {}

with open("igrf13coeffs.txt") as coef_file:
    for i in range(3):
        coef_file.readline()
    line = coef_file.readline()
    years = []
    for year in line.split()[3:-1]:
        years.append(float(year))
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

for degree in h:
    h[degree][0] = [0] * len(years)
    h_sv[degree][0] = 0

print("Years:", years)

print("n m g gsv h hsv")
for degree in g:
    for order in g[degree]:
        print(
            degree,
            order,
            g[degree][order][-1],
            g_sv[degree][order],
            h[degree][order][-1],
            h_sv[degree][order],
        )
