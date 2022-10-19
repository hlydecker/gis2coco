# utils.py

import datetime
import json
import os
import re
import fnmatch
from PIL import Image
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio as rio
from rasterio import plot as rasterplot
from osgeo import gdal
import geoplot.crs as gcrs
import matplotlib.pyplot as plt
from shapely.geometry import box


def plot_raster_and_polygons(raster, polygons, export_plot = False, raster_band = 1):

    """
    Plot a raster and polygon later together
    Arguements:
        raster: raster layer to plot
        polygons: polygons that you want to plot on top of the raster
        export_plot: boolean if set to true will save the plot. Defaults to False
        raster_band: integer specifying the raster band to plot
    Returns:
        A plot of with the raster and polygon overlay.
    """
    
    f, ax = plt.subplots()

    geotiff = rio.open(raster)
    # plot raster 
    rasterplot.show(
        geotiff.read(raster_band),  # use tiff.read(1) with your data
        extent=[geotiff.bounds[0], geotiff.bounds[2], geotiff.bounds[1], geotiff.bounds[3]],
        ax=ax,
    )

    # plot shapefiles
    polygons.plot(ax=ax, facecolor='w', edgecolor='b', alpha=0.5)

    if export_plot == True:
        plt.savefig(f"{geotiff}_{polygons}.jpg")
        print()
    plt.show()