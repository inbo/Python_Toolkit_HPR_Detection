#Standard libraries
import math
import os
#Third party libraries
import rasterio
from rasterio.transform import xy
from rasterio import plot as rasterio_plot
from rasterio import warp as rasterio_warp  #calculate_default_transform, reproject, Resampling
import matplotlib.pyplot as plt
import cv2
import numpy as np
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
import geopandas
import shapely
from shapely import affinity
from shapely import Point, LineString
#Local applications

def open_raster_data(filename, target_crs=None):
	raster = rasterio.open(filename)
    if raster.crs != target_crs:
        print(f"Reprojecting raster data to {target_crs}...")
        
        transform, width, height = rasterio_warp.calculate_default_transform(
            raster.crs, target_crs, raster.width, raster.height, *raster.bounds
        )
        kwargs = raster.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        
        with rasterio.open('Example_VITO-reprojected.tif', 'w', **kwargs) as dst:
            for i in range(1, raster.count + 1):
                rasterio_warp.reproject(
                    source=rasterio.band(raster, i),
                    destination=rasterio.band(dst, i),
                    src_transform=raster.transform,
                    src_crs=raster.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=rasterio_warp.Resampling.nearest  # Choose an appropriate resampling
                )
        raster.close()
        raster = rasterio.open('Example_VITO-reprojected.tif')
        print(f"Raster CRS after reprojection: {raster.crs}")