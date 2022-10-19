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


def plot_raster_and_polygons(geotiff, polygons, export_plot = False):

    """
    Plot a raster and polygon later together
    Arguements:
        geotiff: raster layer to plot
        polygons: polygons that you want to plot on top of the raster
        export_plot: boolean if set to true will save the plot. Defaults to False
    Returns:
        A plot of with the raster and polygon overlay.
    """
    
    f, ax = plt.subplots()

    geotiff = rio.open(geotiff)
    # plot raster 
    rasterplot.show(
        geotiff.read(1),  # use tiff.read(1) with your data
        extent=[geotiff.bounds[0], geotiff.bounds[2], geotiff.bounds[1], geotiff.bounds[3]],
        ax=ax,
    )

    # plot shapefiles
    polygons.plot(ax=ax, facecolor='w', edgecolor='b', alpha=0.5)
    #plt.savefig('small_test_4326.jpg')
    print(f"{geotiff.name}_plot.jpg")
    plt.show()
    if export_plot == True:
        plt.savefig(f"{geotiff.name}_polygons.jpg")
