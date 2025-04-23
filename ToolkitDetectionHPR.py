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

Point.__neg__ = lambda self: Point(-self.x, -self.y)
Point.__pos__ = lambda self: self

Point.__add__ = lambda self, other: Point(self.x+other.x, self.y+other.y)
Point.__sub__ = lambda self, other: self.__add__(-other)
Point.__mul__ = lambda self, other: Point(self.x*other, self.y*other)
Point.__truediv__ = lambda self, other: self.__mul__(1/other)

Point.__iadd__ = lambda self, other: self.__add__(other)
Point.__isub__ = lambda self, other: self.__add__(-other)
Point.__imul__ = lambda self, other: self.__mul__(other)
Point.__idiv__ = lambda self, other: self.__truediv__(other)

Point.__eq__ = lambda self, other: ((self.x == other.x) and (self.y == other.y))
Point.__ne__ = lambda self, other: ((self.x != other.x) or (self.y != other.y))

LineString.get_points = lambda self: [Point(coord[0],coord[1]) for coord in self.coords]

class LineEquation():
    def __init__(self, **kwargs):
        if len(kwargs) == 2:
            if 'rho' in kwargs.keys() and 'theta' in kwargs.keys():
                self._a = np.cos(kwargs['theta'])
                self._b = np.sin(kwargs['theta'])
                self._c = -kwargs['rho']
                return

            points = []
            if 'y_intercept' in kwargs.keys():
                points.append(Point(0., kwargs['y_intercept'])) 
            if 'x_intercept' in kwargs.keys():
                points.append(Point(kwargs['x_intercept'], 0.))
            if 'point' in kwargs.keys():
                points.append(kwargs['point'])
            if 'point1' in kwargs.keys():
                points.append(kwargs['point1'])
            if 'point2' in kwargs.keys():
                points.append(kwargs['point2'])


            if 'slope' not in kwargs.keys() and len(points) == 2:
                dp = points[1] - points[0]
                if dp.x != 0.: 
                    kwargs['slope'] = dp.y / dp.x
                else:
                    self._a = 1.
                    self._b = 0.
                    self._c = -points[0].x
                    return

            if 'slope' in kwargs.keys() and len(points) > 0:
                self._a = kwargs['slope']
                self._b = -1.
                self._c = points[0].y - kwargs['slope'] * points[0].x
                return

            raise IOError('The two key-word arguments should be a combination of the cartesian parameters ("slope", "y_intercept", "x_intercept", "point1" and "point2") or the normal representation ("rho" and "theta").')

            

        elif set(kwargs.keys()) == {'a', 'b', 'c'}:
            self._a = kwargs['a']
            self._b = kwargs['b']
            self._c = kwargs['c']

        else:
            raise IOError('The function excepts only two key-word arguments to initialize the LineEquation, except for the key triplet "a", "b", "c".')

    def get_parameters(self):
        return (self._a, self._b, self._c)

    def get_normal_representation(self):
        # one line is determined by multiple line equations
        # first set the line equation to the following form
        # w*cos(theta)*x + w*sin(theta)*y - w*rho = 0 
        # with w a positive real number
        # for rho = 0, theta in the first are fourth quadrant
        """if self._c > 0. or (self._c == 0. and self._a < 0.):
            self._a *= -1.
            self._b *= -1.
            self._c *= -1.
        
        theta = np.arctan2(self._b, self._a)
        rho = -self._c / self._a * np.cos(theta)

        return rho, theta"""
        normal_line = self.get_perpendicular_line(Point(0,0))
        intersection = self.get_intersection_point(normal_line)
        theta = np.arctan2(intersection.y, intersection.x)
        rho = np.linalg.norm([intersection.x, intersection.y])
        return rho, theta

    def get_x_value(self, y):
        return -(self._b*y + self._c) / self._a

    def get_y_value(self, x):
        return -(self._a*x + self._c) / self._b

    def get_perpendicular_line(self, point):
        new_a = self._b
        new_b = -self._a
        new_c = -self._b*point.x + self._a*point.y
        return LineEquation(a=new_a, b=new_b, c=new_c)

    def get_intersection_point(self, other_line):
        a1, b1, c1 = self.get_parameters()
        a2, b2, c2 = other_line.get_parameters()
        if a1/b1 == a2/b2:
            return None
        else:
            x = (c1/b1 - c2/b2) / (-a1/b1 + a2/b2)
            y = -(a1*x + c1) / b1
            return Point(x, y)


class LineSegmentDetector():
    def __init__(self, config={}):
        self.__default_config = {'rho_tolerance': 10.,
                                 'theta_tolerance': np.deg2rad(5.),
                                 'Hough lines': {'rho_resolution': 1,
                                                 'theta_resolution': np.deg2rad(.1),
                                                 'threshold': 30, 
                                                 'minLineLength': 50., 
                                                 'maxLineGap': 30.,
                                                }
                                }
        self._config = self.merge_config(config, self.__default_config)

    def set_binary_raster(self, binary_raster):
        self._binary_raster = binary_raster

    def set_config(self, new_config):
        self._config = merge_config(new_config, self._config)

    def merge_config(self, user_config, default_config):
        """
        Merge the user configuration dictionary with the default configuration dictionary recursively.

        Parameters
        ----------
        user : dict
            The user configuration dictionary.
        default : dict
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
                    user_config[key] = self.merge_config(user_config[key], default_value)
        return user_config

    def run(self, merge_lines=True):
        HLP_config = self._config['Hough lines']
        line_segments = cv2.HoughLinesP(self._binary_raster, HLP_config['rho_resolution'], HLP_config['theta_resolution'],
            threshold=HLP_config['threshold'], minLineLength=HLP_config['minLineLength'], maxLineGap=HLP_config['maxLineGap'])

        # remove unused dimension in numpy array
        line_segments = line_segments[:,0,:]
        # extra dimention to determine point with x and y
        line_segments = np.reshape(line_segments, (len(line_segments),len(line_segments[0])//2,2))
        # convert into list holding LineString objects
        line_segments = [LineString(line_segment) for line_segment in line_segments]

        if merge_lines:
            line_segments = self.merge_similar_line_segments(line_segments, self._config['rho_tolerance'], self._config['theta_tolerance'])

        return line_segments


    def merge_similar_line_segments(self,line_segments, rho_tolerance, theta_tolerance):
        """
        Merges similar line segments into one infinite line
        based on their perpendicular absolute distance (0 < rho)
        and the polar angle theta of the perpendicular line 
        through the origin (-pi < theta < +pi). 
        
        Parameters
        ----------
        line_segments : list of numpy.ndarray
            A list of line segments, where each line_segment 
            is represented as [[x1, y1], [x2, y2]]
        rho_tolerance : float
            Tolerance for the perpendicular distance (rho)
        theta_tolerance : float
            Tolerance for the polar angle (theta) in radians.

        Returns
        -------
        merged_lines : list of numpy.ndarray
            A list of merged similar lines (average rho and theta).
        """

        merged_lines = []
        processed = [False] * len(line_segments)

        for i in range(len(line_segments)):
            if not processed[i]:
                start_point_i, end_point_i = line_segments[i].get_points()
                line_i = LineEquation(point1=start_point_i, point2=end_point_i)
                current_rho, current_theta = line_i.get_normal_representation()
                similar_lines = {'normal':[[current_rho, current_theta]],
                                 'xrange':[start_point_i.x,end_point_i.x]}
                processed[i] = True
                
                for j in range(i + 1, len(line_segments)):
                    if not processed[j]:
                        start_point_j, end_point_j = line_segments[j].get_points()
                        line_j = LineEquation(point1=start_point_j, point2=end_point_j)
                        rho, theta = line_j.get_normal_representation()

                        # check for similar lines
                        if ((abs(theta - current_theta) < theta_tolerance)
                                and (abs(rho - current_rho) < rho_tolerance)):
                            similar_lines['normal'].append([rho, theta])
                            processed[j] = True
                            
                        # check for lines with similar rho and theta
                        # at the -pi and +pi boundary
                        elif (((2*np.pi - abs(theta - current_theta)) < theta_tolerance)
                                and (abs(rho - current_rho) < rho_tolerance)):
                            continue
                            similar_lines['normal'].append([rho, theta + np.sign(current_theta)*2*np.pi])
                            processed[j] = True

                        # check for line with opposite theta
                        elif (((np.pi - abs(theta - current_theta)) < theta_tolerance)
                                and (rho + current_rho < rho_tolerance)):
                            continue
                            similar_lines['normal'].append([-rho, theta + np.sign(current_theta)*np.pi])
                            processed[j] = True

                        if processed[j]:
                            similar_lines['xrange'].extend((start_point_j.x, end_point_j.x))

                # Calculate the average rho and theta
                similar_lines['normal'] = np.array(similar_lines['normal'])
                avg_rho = np.median(similar_lines['normal'][:,0])
                avg_theta = np.median(similar_lines['normal'][:,1])
                x1 = np.min(similar_lines['xrange'])
                x2 = np.max(similar_lines['xrange'])
                y1, y2 = LineEquation(rho=avg_rho, theta=avg_theta).get_y_value(np.array([x1, x2]))
                merged_lines.append(LineString([[x1, y1], [x2, y2]]))

        return merged_lines

    