#Standard libraries
from typing import Optional
#Third party libraries
import numpy as np
#Local applications

class DitchesDetector():
    def __init__(self):
        pass

    def has_plot_ditches(self, plot, relief):
        pass

    def select_image_background(
        self, image, 
        background_estimation_method : str='median', 
        threshold_factor : float=2., 
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
        threshold : float, default=2.
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