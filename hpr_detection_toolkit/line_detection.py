#Standard libraries
import os
#Third party libraries
import cv2
import yaml
import numpy as np
from shapely import Point, LineString
#Local applications
from . import utils

# Extend shapely.Point class to include regular arithmic operations
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
# Extend shapely.Point class to evaluate if points are the same
Point.__eq__ = lambda self, other: ((self.x == other.x) and (self.y == other.y))
Point.__ne__ = lambda self, other: ((self.x != other.x) or (self.y != other.y))
# Extend shapely.LineString class to return the begin and end shapely.Point
LineString.get_points = lambda self: [Point(coord[0],coord[1]) for coord in self.coords]

class LineEquation():
    """
    Class to represent a infinite line in two-dimensional space by the following equation
    a*x + b*y + c = 0

    Attributes
    ----------
    a : float
    b : float
    c : float

    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : dict
            Values that define the 2D line described by key words.
            The following combination are allowed:
            - a: float, b: float, c: float (general form)
            - rho: float, theta: float (normal form)
            - x_intercept: float, slope: float
            - y_intercept: float, slope: float
            - x_intercept: float, y_intercept: float
            - point: shapely.Point, slope: float
            - point1: shapely.Point, point2: shapely.Point
        """

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
        """
        return the parameters a, b and c of the general line equation

        Returns
        -------
        (a, b, c) : tuple of 3 floats
        """

        return (self._a, self._b, self._c)

    def get_normal_representation(self):
        """
        return the parameters rho and theta of the line's normal representation

        Returns
        -------
        (rho, theta) : tuple of 2 floats
        """

        normal_line = self.get_perpendicular_line(Point(0,0))
        intersection = self.get_intersection_point(normal_line)
        theta = np.arctan2(intersection.y, intersection.x)
        rho = np.linalg.norm([intersection.x, intersection.y])
        return rho, theta

    def get_x_value(self, y):
        """
        evaluate the line equation for x based on y-value(s)

        Returns
        -------
        x : float or 1D numpy array
        """
        return -(self._b*y + self._c) / self._a

    def get_y_value(self, x):
        """
        evaluate the line equation for y based on x-value(s)

        Returns
        -------
        y : float or 1D numpy array
        """
        return -(self._a*x + self._c) / self._b

    def get_perpendicular_line(self, point):
        """
        calculate the line perpendicular to this line and going through a certain point

        Parameters
        ----------
        point : shapely.Point
            point through which the perpendicular line has to go

        Returns
        -------
        perp_line : LineEquation
        """

        new_a = self._b
        new_b = -self._a
        new_c = -self._b*point.x + self._a*point.y
        return LineEquation(a=new_a, b=new_b, c=new_c)

    def get_intersection_point(self, other_line):
        """
        calculate the intersection point of this line with another line

        Parameters
        ----------
        other_line : LineEquation

        Returns
        -------
        intersection : shapely.Point
        """
        
        a1, b1, c1 = self.get_parameters()
        a2, b2, c2 = other_line.get_parameters()
        if a1/b1 == a2/b2:  # check if lines are parallel
            return None
        else:
            x = (c1/b1 - c2/b2) / (-a1/b1 + a2/b2)
            y = -(a1*x + c1) / b1
            return Point(x, y)


class LineSegmentDetector():
    """
    Detects lines in binary 2D images
    """

    def __init__(self, config={}):
        """
        Parameters
        ----------
        config : dict, default={}
            User configuration for the line segment detection.
            When empty, the default configuration is used.
            Otherwise the default configuration is overwritten
            for key words present in the config file.
        """

        module_directory = os.path.dirname(os.path.abspath(__file__))
        config_filename = os.path.join(module_directory, 'config_default.yaml')
        self.__default_config = yaml.safe_load(open(config_filename,'r'))['line_segment_detection']
        self._config = utils.merge_config(config, self.__default_config)

    def process(self, binary_image, merge_lines=False):
        """
        process a binary image and extract line segments

        Parameters
        ----------
        binary_image : 2D numpy.array with dtype = numpy.uint8
            raster image of zeros and ones.
        merge_lines : boolean, default=False
            merge smaller line segments with similar line equationa
            into one line segments.

        Returns
        -------
        line_segments : list of shapely.LineString
            detected line segments in vector format 
        """

        self._line_segments = None
        
        if binary_image.dtype != np.uint8:
            raise TypeError('The dtype of the binary image should be numpy.uint8')

        HLP_config = self._config['Hough lines']
        line_segments = cv2.HoughLinesP(binary_image, HLP_config['rho_resolution'], np.deg2rad(HLP_config['theta_resolution']),
            threshold=HLP_config['threshold'], minLineLength=HLP_config['minLineLength'], maxLineGap=HLP_config['maxLineGap'])
        
        if line_segments is not None:
            line_segments = line_segments[:,0,:]  # remove unused dimension in numpy array
            line_segments = np.reshape(line_segments, (len(line_segments),len(line_segments[0])//2,2))  # extra dimention to represent point as [x, y]
            line_segments = [LineString(line_segment) for line_segment in line_segments]  # convert into list holding LineString objects

            if merge_lines:
                line_segments = self._merge_similar_line_segments(line_segments, self._config['rho_tolerance'], np.deg2rad(self._config['theta_tolerance']))

            self._line_segments = line_segments

    def get_line_segments(self):
        """
        return the line segments detected in the binary image

        Returns
        -------
        line_segment : list of shapely.LineString
        """
        return self._line_segments

    def _merge_similar_line_segments(self,line_segments, rho_tolerance, theta_tolerance):
        """
        Merges smaller similar line segments into bigger segments
        based on their perpendicular absolute distance (0 < rho)
        and the polar angle theta of the perpendicular line 
        through the origin (-pi < theta < +pi). 
        
        Parameters
        ----------
        line_segments : list of shapely.LineString
            List of detected line segments.
        rho_tolerance : float
            Tolerance for the perpendicular distance (rho)
            to evaluate similarity between line segments.
        theta_tolerance : float
            Tolerance for the polar angle (theta) in radians
            to evaluate similarity between line segments.

        Returns
        -------
        merged_lines : list of shapely.LineString
            A list of the merged line segments.
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