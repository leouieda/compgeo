"""
Week 7: Make grids of IGRF
"""
import datetime
import numpy as np
import xarray as xr
import verde as vd
import pygmt
import igrf


def igrf_grid(region, spacing, height, date):
    """
    Make a grid of Bx, By, Bz at a uniform height.
    """
    longitude, latitude = vd.grid_coordinates(region, spacing=spacing, meshgrid=False)

    shape = (latitude.size, longitude.size)
    be = np.zeros(shape)
    bn = np.zeros(shape)
    bu = np.zeros(shape)

    for i, lat in enumerate(latitude):
        for j, lon in enumerate(longitude):
            be_point, bn_point, bu_point = igrf.igrf(lon, lat, height, date)
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


if __name__ == "__main__":
    grid = igrf_grid((0, 360, -90, 90), spacing=3, height=1000, date=datetime.datetime.today())
    print(grid)
    grid["amplitude"] = np.sqrt(grid.be**2 + grid.bn**2 + grid.bu**2)
    fig = pygmt.Figure()
    fig.grdimage(grid.amplitude, projection="W0/20c", cmap="viridis", frame=True)
    fig.colorbar()
    fig.coast(shorelines=True)
    fig.show()
