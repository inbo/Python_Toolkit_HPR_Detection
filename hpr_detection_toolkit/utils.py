#Standard libraries
import math
import os
#Third party libraries
import rasterio
import rasterio.warp
import rasterio.vrt
import cv2
import numpy as np
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
import geopandas
import shapely
from shapely.affinity import affine_transform
from shapely import Point, LineString
#Local applications


def open_raster_data(filename, target_crs=None):
    try:
        raster = rasterio.open(filename)
    except rasterio.RasterioIOError:
        raise IOError("Error: Could not open or read the TIF image.")

    if target_crs is not None and raster.crs != target_crs:
        print(f"Reprojecting raster data to {target_crs}...")
        
        transform, width, height = rasterio.warp.calculate_default_transform(
            raster.crs, target_crs, raster.width, raster.height, *raster.bounds
        )
        kwargs = raster.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        
        filename_base, filename_ext = os.path.splitext(filename)
        with rasterio.open(f'{filename_base}-reprojected{filename_ext}', 'w', **kwargs) as output_file:
            for i in range(1, raster.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(raster, i),
                    destination=rasterio.band(output_file, i),
                    src_transform=raster.transform,
                    src_crs=raster.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=rasterio.warp.Resampling.nearest  # Choose an appropriate resampling
                )
        raster.close()
        raster = rasterio.open(f'{filename_base}-reprojected{filename_ext}')
        print(f"Raster CRS after reprojection: {raster.crs}")

    else:
        print("Raster data in original CRS.")

    return raster

def print_raster_metadata(raster):
    print("-"*40)
    print(f"METADATA - {raster.name}")
    print("-"*40)
    print(f"Image bands:\t{raster.count}")
    print(f"Image shape:\t({raster.width},{raster.height})")
    print(f"CRS:\t\t{raster.crs}")
    print(f"Bounds:\t\tleft-right ({raster.bounds.left}, {raster.bounds.right})")
    print(f"\t\ttop-bottom ({raster.bounds.top}, {raster.bounds.bottom})")
    print("")
    print(f"Geotransform:\n{raster.transform}")
    print("-"*40)

def pixel_to_georef(geometry, transform_matrix):
    try:
        return [affine_transform(geom, transform_matrix) for geom in geometry]
    except TypeError:
        pass
    else:
        return affine_transform(geometry, transform_matrix)