#Standard libraries
import os
from typing import Optional
#Third party libraries
import yaml
import numpy as np
import geopandas
import shapely
from scipy.ndimage import gaussian_filter
#Local applications
from .line_detection import LineSegmentDetector
from . import utils


class HprDitchDetector():
    """ Detect linear structures in a raster map 
    based on microrelief within a plot of land """
    def __init__(self, relief, config={}):
        """
        Attributes
        ----------
        relief : rasterio.DataReader-like object
            The raster image in gray-scale where linear structures in high or low values
            represent ditches in the landplot. Can be microrelief or processed relief map
            such as the map made by VITO.
        config : dict, default={}
            User's configuration for the hpr ditch detection.
            When empty, the default configuration is used.
            Otherwise the default configuration is overwritten
            for key words present in the user's configuration.
        """
        self._LSD = LineSegmentDetector()
        self._relief = relief
        self._crs = relief.crs

        module_directory = os.path.dirname(os.path.abspath(__file__))
        config_filename = os.path.join(module_directory, 'config_default.yaml')
        self.__default_config = yaml.safe_load(open(config_filename,'r'))['hpr_ditch_detection']
        self._config = utils.merge_config(config, self.__default_config)

    def process(self, landplot):
        """
        Detect ditches and determine the percentage of hpr covered land

        Parameters
        ----------
        landplot : geopandas.DataFrame
            Georeferenced list holding the landplot
        """

        # Make sure results of previous processing runs are erased
        self._ditches = None
        self._buffer_zone = None
        self._hpr_fraction = None

        # Perform necessary processes
        self._ditches = self._detect_ditches_in_plot(landplot, relief=self._relief, merge_lines=self._config['merge_lines'])
        self._buffer_zone = self._buffer_ditches(self._ditches, landplot, **self._config['buffer_zone'])
        self._hpr_fraction = self._buffer_zone.iloc[0].geometry.area / landplot.iloc[0].geometry.area

    def get_hpr_fraction(self):
        """
        return the percentage of landplot surface that is covered by the buffer zone

        Returns
        -------
        hpr_fraction : float
        """
        return self._hpr_fraction

    def get_ditches(self, multilinestring=False):
        """
        return the location of ditches represented by line segments

        Returns
        -------
        ditches : list of shapely.LineString
        """
        if not multilinestring:
            return self._ditches
        else:
            return geopandas.GeoDataFrame(geometry=[shapely.MultiLineString(list(self._ditches.geometry))], crs=self._crs)

    def get_buffer_zone(self):
        """
        return the buffer zone around the ditches

        Returns
        -------
        buffer_zone : geopandas.DataFrame
        """
        return self._buffer_zone

    def _detect_ditches_in_plot(self, landplot, relief, merge_lines=False):
        """
        detected ditches (linear structures) in the relief raster image 

        Parameters
        ----------
        landplot : geopandas.DataFrame
            Georeferenced list holding the landplot
        relief : shapely.DataReader-like object
            The raster image in gray-scale where linear structures in high or low values
            represent ditches in the landplot. Can be microrelief or processed relief map
            such as the map made by VITO.
        merge_lines : boolean, default=False
            If true, line segments with a similar line equation and in the vicinity of
            each other or merged into one single line segment.

        Returns
        -------
        geopandas.DataFrame
            Georeferenced list of all the detected line segments
        """

        clipped_relief = utils.clip_raster(relief, landplot.iloc[0].geometry.geoms)  # only use pixels within landplot

        image = clipped_relief.read(1)
        mask_background = self._select_image_background(image, **self._config['filter_background'])
        image_masked = np.where(mask_background, 0., image)
        image_binary = np.where(mask_background, image_masked, 1.)

        self._LSD.process(image_binary.astype(np.uint8), merge_lines=merge_lines)
        line_segments = self._LSD.get_line_segments()

        line_segments_geometry = utils.pixel_to_georef(line_segments, clipped_relief.transform.to_shapely())
        return geopandas.GeoDataFrame(geometry=line_segments_geometry, crs=self._crs)

    def _buffer_ditches(self, lines, landplot, distance):
        """
        create a buffer polygon around the ditches within a landplot

        Parameters
        ----------
        lines : geopandas.DataFrame
            Georeferenced list of line segments representing ditches within the landplot
        landplot : geopandas.DataFrame
            Georeferenced list holding the landplot
        distance : float
            Buffer distance around the line in units of crs

        Returns
        -------
        geopandas.DataFrame
            Georeferenced list holding the buffer zone within the landplot
        """

        buffer_gdf = lines.copy()  # Create a copy to store the buffer geometries
        buffer_gdf['geometry'] = lines.geometry.buffer(distance) # directly assign
        buffer_union_gdf = geopandas.GeoDataFrame(geometry=[buffer_gdf.union_all()], crs=self._crs)        
        return geopandas.overlay(buffer_union_gdf, landplot, how='intersection')

    def _select_image_background(
        self, image, 
        background_estimation_method : str='median', 
        threshold_factor : float=3., 
        include_high : bool=True, 
        include_low : bool=True, 
        return_background_estimate : bool=False,
        gaussian_sigma : Optional[float]=None
    ):
        """
        Estimate the background value for every pixel in the image.


        Parameters
        ----------
        image : 2D numpy.ndarray 
            The input image (grayscale).
        background_estimation_method : str, default='mean' 
            The method used to estimate the background. Options are: 'mean', 'median', 'gaussian'.
        threshold : float, default=3.
            The number of standard deviations by which a pixel's value must differ from the 
            background to be considered a high-value or low-value pixel.
        gaussian_sigma : float, optional, default=None: 
            The standard deviation for the Gaussian filter if 'gaussian' is chosen as the 
            background estimation method. If None, a reasonable value based on image size 
            might be used.


        Returns
        -------
        mask_background : numpy.ndarray
            A boolean mask of the same shape as the input image, 
            where True indicates a background pixel and False otherwise.
        background_estimate : numpy.ndarray, optional
            An image of the same shape as the input image holding
            the estimated background value for each pixel.
        """

        image = np.asarray(image)

        # Grayscale image
        if background_estimation_method == 'mean':
            background = np.ones_like(image)*np.mean(image)
        elif background_estimation_method == 'median':
            background = np.ones_like(image)*np.median(image)
        elif background_estimation_method == 'gaussian':
            if gaussian_sigma is None:
                gaussian_sigma = max(image.shape) // 50  # Heuristic for sigma
            background = gaussian_filter(image, sigma=gaussian_sigma)
        else:
            raise ValueError(f"Invalid background_estimation_method: {background_estimation_method}")

        standard_deviation = np.std(image - background, mean=0.)
        high_value_mask = (image > (background + threshold_factor * standard_deviation))
        low_value_mask = (image < (background - threshold_factor * standard_deviation))
        
        mask_background = ~((high_value_mask & include_high) | (low_value_mask & include_low))

        if return_background_estimate:
            return mask_background, background
        else:
            return mask_background

    def adjust_config(self, adjust_config):
        """
        Adjust the set configuration.

        Parameters
        ----------
        adjust_config : dict
            Configuration that should be adjusted.

        """
        self._config = utils.merge_config(adjust_config, self._config)

    def reset_config(self, new_config={}):
        """
        Reset the configuration to default except for the configuration given

        Parameters
        ----------
        new_config : dict
            Configuration that should be set differently from the default.

        """
        self._config = utils.merge_config(new_config, self._default_config)