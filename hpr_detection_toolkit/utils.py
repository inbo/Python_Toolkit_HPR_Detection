#Standard libraries
#Third party libraries
import rasterio
import rasterio.warp
import rasterio.vrt
import rasterio.mask
import rasterio.io
import numpy as np
import geopandas
import shapely
import shapely.affinity
from shapely import Point, LineString
#Local applications


def open_raster_data(filename, target_crs=None):
    try:
        raster = rasterio.open(filename)
    except rasterio.RasterioIOError:
        raise IOError("Error: Could not open or read the TIF image.")

    if target_crs is not None and raster.crs != target_crs:
        print(f"Reprojecting raster data to {target_crs}...")
        raster = rasterio.vrt.WarpedVRT(raster, crs=target_crs)
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
        return [shapely.affinity.affine_transform(geom, transform_matrix) for geom in geometry]
    except TypeError:
        pass
    else:
        return shapely.affinity.affine_transform(geometry, transform_matrix)


def merge_config(user_config, default_config):
    """
    Merge the user configuration dictionary with the default configuration dictionary recursively.

    Parameters
    ----------
    user_config : dict
        The user configuration dictionary.
    default_config : dict
        The default configuration dictionary.

    Returns
    -------
    dict
        The merged configuration dictionary.

    """
    if isinstance(user_config, dict) and isinstance(default_config, dict):
        for key, default_value in default_config.items():
            if key not in user_config:
                user_config[key] = default_value
            else:
                user_config[key] = merge_config(user_config[key], default_value)
    return user_config


def clip_raster(raster, geoms):
    """
    Clip a raster image within a rasterio.DataReader object

    Parameters
    ----------
    raster : rasterio.DataReader-like object
        The original raster image including the metadata.
    geoms : shapely geometry sequence
        List holding the polygon(s) used as clipping mask

    Returns
    -------
    clipped_raster : rasterio.DataReader-like object
        The clipped raster image including the metadata.

    """

    # Select only the pixels from the relief map within the landplot
    out_image, out_transform = rasterio.mask.mask(raster, geoms, crop=True)
    out_meta = raster.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform,
    })

    # Make a rasterio.DataReader-like object to manipulate the clipped image
    # without saving the clipped image to a temporary file and opening it
    memfile = rasterio.io.MemoryFile()
    with memfile.open(**out_meta) as dataset:  # Open memory-file in write modus
            dataset.write(out_image)  # Write the NumPy array to the dataset
    clipped_raster = memfile.open()  # Open memory-file in read modus

    return clipped_raster