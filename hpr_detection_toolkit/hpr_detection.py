#Standard libraries
from typing import Optional
#Third party libraries
import numpy as np
import geopandas
#Local applications
from Python_Toolkit_HPR_Detection.hpr_detection_toolkit.line_detection import LineSegmentDetector
from Python_Toolkit_HPR_Detection.hpr_detection_toolkit import utils

class HprDitchDetector():
    def __init__(self):
        self._LSD = LineSegmentDetector()
        self._plot = None
        self._relief = None
        self._crs = None
        self._ditches = None
        self._buffer_zone = None
        self._hpr_fraction = None

    def begin(self, plot, relief):
        self._plot = plot
        self._relief = relief
        self._crs = relief.crs

    def run(self, threshold_factor=3., merge_lines=False):
        self._ditches = self._detect_ditches_in_plot(threshold_factor=threshold_factor, merge_lines=merge_lines)
        self._buffer_zone = self._buffer_ditches(self._ditches)
        self._hpr_fraction = self._buffer_zone.iloc[0].geometry.area / self._plot.iloc[0].geometry.area

    def end(self):
        self._plot = None
        self._relief = None
        self._crs = None
        self._ditches = None
        self._buffer_zone = None
        self._hpr_fraction = None

    def get_hpr_fraction(self):
        return self._hpr_fraction

    def get_ditches(self):
        return self._ditches

    def get_buffer_zone(self):
        return self._buffer_zone

    def _detect_ditches_in_plot(self, relief=None, threshold_factor=3., merge_lines=False):
        if relief is None:
            relief = self._relief

        image = relief.read(1)
        mask_background = self._select_image_background(image, threshold_factor=threshold_factor)
        image_masked = np.where(mask_background, 0., image)
        image_binary = np.where(mask_background, image_masked, 1.)

        self._LSD.set_binary_raster(image_binary.astype(np.uint8))
        self._LSD.run(merge_lines=merge_lines)
        line_segments = self._LSD.get_line_segments()

        line_segments_geometry = utils.pixel_to_georef(line_segments, relief.transform.to_shapely())
        return geopandas.GeoDataFrame(geometry=line_segments_geometry, crs=self._crs)

    def _buffer_ditches(self, lines, plot=None, buffer_distance=30.):
        if plot is None:
            plot = self._plot

        buffer_gdf = lines.copy()  # Create a copy to store the buffer geometries
        buffer_gdf['geometry'] = lines.geometry.buffer(buffer_distance) # directly assign
        buffer_union_gdf = geopandas.GeoDataFrame(geometry=[buffer_gdf.union_all()], crs=self._crs)        
        return geopandas.overlay(buffer_union_gdf, plot, how='intersection')

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
        estimated background value for every pixel in the image.

        Parameters
        ----------
        image : numpy.ndarray 
            The input image (grayscale).
        background_estimation_method : str, default='mean': 
            The method used to estimate the background. Options are: 'mean', 'median', 'gaussian'.
        threshold : float, default=3.
            The number of standard deviations by which a pixel's value must differ from the 
            background to be considered a high-value pixel.
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