# Copyright (C) 2021-2023 Geoscience Australia
# 
# The minpot-toolkit is released under the Apache License, Version 2.0 
# (the "License");you may not use this software except in compliance with 
# the License. You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# The project uses third party components which may have different licenses. 
# Please refer to individual components for more details.
"""
Feature Engineering module
-----------------
This module is used to create mappable proxies from geophysical models.
    
"""

import numpy as np
#import numba, time
import time
import toolkit.functions as fn
import toolkit.plotting as pt
from scipy.interpolate import RegularGridInterpolator
from scipy import signal

class FeatureManager():    
    def __init__(self, commandList, features):
        """
        A class to facilitate the generation of features from input geophysical 
        models, given input commands from xml configuration file.
        
        
        Parameters
        ----------
        model_dict: dictionary
            A dictionary of the geophysical models which are available for use.
        
        commandList: list
            A list of commands read from XML config file.
        
        features: xml.ElementTree node
            Nodes of the XML config file which describe features that are to be
            generated.
            
        
        Calls
        -----
        functions.changetype
        
        
        """
        self.features = features
        
        for row in commandList:
            if hasattr(self, row[0]):
                tp = type(getattr(self, row[0]))
                val = fn.changetype(row[1], tp)
                setattr(self, row[0], val)
            #end if
        #end for
        
        #Compile distance calculation function
        a = np.array([0, 1]).astype("float32")
        _ = cosdist_to_feat(a, a, a, a, a, a, a, a, a, a)
        _ = sindist_to_feat(a, a, a, a, a, a, a, a, a, a, a, a)
    #end func
    
    def Process(self, model_dict):
        """
        Generates a "feature" object for each feature defined in the xml config 
        file.
        
        
        Calls
        -----
        feature_object
        
        feature_object.operate_on_existing_features
        
        feature_object.generate_feature
        
        feature_object.distance_to_feature
        
        
        """
        self.feature_dict = {}
        for feature in self.features:
            t0 = time.time()
            key = feature.attrib["name"]
            if "use_existing_features" in feature.attrib:
                use_existing_features = \
                    fn.changetype(feature.attrib["use_existing_features"], 
                                  bool)
            else:
                use_existing_features = False
            #end if
            if use_existing_features:
                feat = feature_object(feature)
                feat.operate_on_existing_features(self.feature_dict)
                feat.distance_to_feature()
            else:
                feat = feature_object(feature)
                feat.generate_feature(model_dict)
                feat.distance_to_feature()
            #end if
            self.feature_dict[key] = feat
            print(key, time.time() - t0)
        #end for
    #end func    
#end class
    
class feature_object():    
    def __init__(self, feature):
        """
        A class to facilitate the engineering of a feature using input 
        geophysical models. 
        
        
        Parameters
        ----------
        feature: xml.ElementTree node
            A node which contains information and criteria required to generate
            the feature from the input geophysical models.
            
        model_dict: dictionary
            A dictionary containing the required input geophysical models. 
            Optional if feature is to be generated from already existing 
            features.
            
        feature_dict: dictionary
            A dictionary containing already generated features. Optional if the
            feature has criteria which only depend on geophysical models.
            
        generate_plots: boolean
            Flag to allow the generation of a figure which represents the 
            feature.
            
            
        Calls
        -----
        functions.changetype
        
        
        """
        self.feature = feature
        self.name = feature.attrib["name"]
        if "generate_plots" in feature.attrib:
            self.generate_plots = \
                fn.changetype(feature.attrib["generate_plots"], bool)
        else:
            self.generate_plots = False
        #end if
    #end func
    
    def generate_feature(self, model_dict):
        """
        Routine to calculate function values for the type of feature under 
        investigation e.g. original values, gradient, etc. 
        
        
        Calls
        -----
        functions.changetype
        
        project_onto_surface
        
        project_from_surface
        
        truth_array_from_gradient_sph
        
        truth_array_from_values
        
        
        """
        
        arrays = []
        isfinite_list = []
        criteria = self.feature.findall("criterion")
        truth_list = list()
        for criterion in criteria:
            truth_list.append("project_to_surface" in criterion.attrib)
        #end for
        if np.all(truth_list) or ~np.any(truth_list):
            proj_from_surface = False
        else:
            proj_from_surface = True
        #end if
        for criterion in criteria:
            name = criterion.attrib["name"]
            model_type = criterion.attrib["model_type"]
            attribute, target_range, target_value = \
                fn.changetype(criterion.attrib["description"], list, dtype=str)
            proj_to_surface = "project_to_surface" in criterion.attrib
            
            model = model_dict[model_type]
            print(name, model_type, attribute, target_range, target_value)
            self.gclon, self.gclat = np.meshgrid(model['longitude'], 
                                                 model['latitude'])
            self.gcz = model['depth']
            values = model['values']
            
            if attribute == 'Value':
                arr, isfinite = truth_array_from_values(values, target_range, 
                                                        target_value)
            elif attribute == 'Gradient':
                arr, isfinite = \
                    truth_array_from_gradient_sph(self.gclon, self.gclat, 
                                                  self.gcz, values, 
                                                  target_range, target_value)
            #end if
            
            if proj_to_surface:
                min_depth, max_depth = \
                    fn.changetype(criterion.attrib["project_to_surface"], 
                                  "array", dtype=float)
                arr = project_onto_surface(arr, self.gcz, min_depth=min_depth,
                                           max_depth=max_depth)
                isfinite = project_onto_surface(isfinite, self.gcz, 
                                                min_depth=min_depth, 
                                                max_depth=max_depth)
            #end if
            if proj_from_surface:
                arr = project_from_surface(arr, self.gcz)
                isfinite = project_from_surface(isfinite, self.gcz)
            #end if
            
            arrays.append(arr)
            isfinite_list.append(isfinite)
        #end if
        vals = np.all(arrays, axis=0)*1.0
        if vals.shape[-1] == 1:
            self.gcz = np.array([0])
        #end if
        isfinite = np.all(isfinite_list, axis=0)
        vals[~isfinite] = np.nan
        self.feature_values = vals
    #end func
    
    def operate_on_existing_features(self, feature_dict):
        """
        Routine to examine existing features and combine characteristics of 
        them, e.g. to look at magnetite/haematite boundaries. 
        
        
        Calls
        -----
        project_onto_surface
        
        project_from_surface
        
        truth_array_from_values
        
        
        """
        
        arrays = []
        isfinite_list = []
        criteria = self.feature.findall("criterion")
        truth_list = list()
        for criterion in criteria:
            truth_list.append("project_to_surface" in criterion.attrib)
        #end for
        if np.all(truth_list) or ~np.any(truth_list):
            proj_from_surface = False
        else:
            proj_from_surface = True
        #end if
        for criterion in criteria:
            name = criterion.attrib["name"]
            feature_name = criterion.attrib["feature_name"]
            value, target_range, target_value = \
                fn.changetype(criterion.attrib["description"], list, dtype=str)
            proj_to_surface = "project_to_surface" in criterion.attrib
            
            model = feature_dict[feature_name]
            if not hasattr(self, 'gcz'):
                self.gclon, self.gclat = model.gclon, model.gclat
                self.gcz = model.gcz
            #end if
            
            print(name, feature_name, 'distance', target_range, target_value)
            
            distances = model.distances
            
            arr, isfinite = \
                truth_array_from_values(distances, target_range, target_value)
            
            if proj_to_surface:
                min_depth, max_depth = \
                    fn.changetype(criterion.attrib["project_to_surface"], 
                                  "array", dtype=float)
                arr = project_onto_surface(arr, self.gcz, min_depth=min_depth,
                                           max_depth=max_depth)
                isfinite = project_onto_surface(isfinite, self.gcz, 
                                                min_depth=min_depth, 
                                                max_depth=max_depth)
            #end if
            if proj_from_surface:
                arr = project_from_surface(arr, self.gcz)
                isfinite = project_from_surface(isfinite, self.gcz)
            #end if
            
            arrays.append(arr)
            isfinite_list.append(isfinite)
        #end if
        vals = np.all(arrays, axis=0)*1.0
        if vals.shape[-1] == 1:
            self.gcz = np.array([0])
        #end if
        isfinite = np.all(isfinite_list, axis=0)
        vals[~isfinite] = np.nan
        self.feature_values = vals
    #end func
    
    def distance_to_feature(self):
        """
        Finds the points on the edges of the feature, and computes the great
        circle distance of every grid point to the closest point on the edge of
        the feature.
        
        
        Calls
        -----
        gcdist_to_feat
        
        functions.gradient
        
        
        """
        lon = self.gclon
        lat = self.gclat
        if (lon[0,1] - lon[0,0] <= 0.025) or (lat[1,0] - lat[0,0] <= 0.025):
            spacing = "small"
        else:
            spacing = "large"
        #end if
        depth = self.gcz
        vals = self.feature_values
        isnan = np.isnan(vals)
        arr = (vals > 0)
        grad = fn.gradient(lon, lat, depth, vals)
        grad = (grad > 0)
        grad = ~arr*grad
        kernel = np.array([[[0], [1], [0]], 
                           [[1], [1], [1]], 
                           [[0], [1], [0]]])
        grad = signal.fftconvolve(grad*1, kernel, mode="same")
        grad = (grad >= 0.5)
        grad = arr*grad
        distances = np.empty((lon.shape[0], lon.shape[1], depth.shape[0]))
        for i in range(len(depth)):
            gradi = grad[:,:,i]
            edges = np.transpose(np.array([lon[gradi], lat[gradi]]))
            points = np.vstack([lon.flatten(), lat.flatten()]).T
            dist = gcdist_to_feat(points, edges, spacing=spacing)
            dist = np.reshape(dist, gradi.shape)
            dist[arr[:,:,i]] = -dist[arr[:,:,i]]
            distances[:,:,i] = dist
        #end for
        distances = np.transpose(distances, axes=(1, 0, 2))
        x = lon[0,:]
        y = lat[:,0]
        if len(depth) == 1:
            fun = RegularGridInterpolator((x, y), distances[:,:,0], 
                                          bounds_error=False, fill_value=None)
        else:
            fun = RegularGridInterpolator((x, y, depth), distances,
                                          bounds_error=False, fill_value=None)
        #end if
        self.fun = fun
        self.distances = np.copy(np.transpose(distances, axes=(1, 0, 2)))
        self.distances[isnan] = np.nan
    #end func
    
    def Plots(self, points_dict, wd_images, bgimage=None):
        """
        Function to generate figure representing feature.
        
        
        Parameters
        ----------
        points_dict: dictionary
            Dictionary which must contain 'point_data_filt', an array of 
            deposit locations.
            
        wd_images: string
            Path in which to save the figure.
            
        bgimage: string, or matplotlib.imread object
            Background image to overlay the feature figure onto.
            
            
        Calls
        -----
        plotting.plot_feature_at_depth_slice
        
        
        """
        if self.generate_plots:
            pt.plot_feature_at_depth_slice(self, points_dict, self.gcz, 
                                           wd_images, bgimage=bgimage)
        #end if
    #end func
    
    def print_description(self):
        """
        Function to print information about the feature to the console, or to
        stdout. This function is called to print feature information to an 
        output file.
        
        
        Calls
        -----
        XML
        
        
        Called by
        ---------
        IO.OutputManager.print_to_file
        
        functions.changetype
        
        
        """
        print('Feature =', self.name)
        print('Model types =', self.feature.attrib["model_types"])
        print('Criteria:')
        criteria = self.feature.findall("criterion")
        for criterion in criteria:
            name = criterion.attrib["name"]
            value, target_range, target_value = \
                fn.changetype(criterion.attrib["description"], list, dtype=str)
            if "model_type" in criterion.attrib:
                model = criterion.attrib["model_type"]
                print('%s:'%name, model, value, target_range, target_value)
            elif "feature_name" in criterion.attrib:
                feature = criterion.attrib["feature_name"]
                print('%s:'%name, value, target_range, target_value, 'from', 
                      feature)
            #end if
            if "project_to_surface" in criterion.attrib:
                low, high = \
                    fn.changetype(criterion.attrib["project_to_surface"], 
                                  list, dtype=str)
                print('Depth slices from', low, 'to', high, 
                      'projected to surface')
            #end if
        #end for
    #end func
#end class
    
def truth_array_from_values(values, target_range, target_value):
    """
    Checks how the values in an array compare to a target value.
    
    
    Parameters
    ----------
    values : numpy array
        3D array of values from a geophysical model
        
    target_range : string
        Either 'less than' or 'greater than'
        
    target_value : float
        Target value to compare values array to.
        
    
    Returns
    -------
    arr : numpy array
        3D array of boolean values representing whether 'values' fulfilled the
        condition applied.
        
    isfinite : numpy array
        3D array of boolean values representing whether or not 'values' 
        contained finite values.
        
    
    """
    if target_range == 'less than':
        arr = values <= float(target_value)
    elif target_range == 'greater than':
        arr = values >= float(target_value)
    #end if
    isfinite = np.isfinite(values)
    return arr, isfinite
#end func
    
def truth_array_from_gradient_sph(lat, lon, z, values, target_range, 
                                  target_value):
    """
    Checks how the the gradient of values in an array compare to a target 
    value.
    
    
    Parameters
    ----------
    values : numpy array
        3D array of values from a geophysical model
        
    target_range : string
        Either 'less than' or 'greater than'
        
    target_value : float
        Target value to compare values array to.
        
    
    Returns
    -------
    arr : numpy array
        3D array of boolean values representing whether 'values' fulfilled the
        condition applied.
        
    isfinite : numpy array
        3D array of boolean values representing whether or not 'values' 
        contained finite values.
        
    
    """
    r = 6371000 - z
    grad = fn.gradient_sph(lon, lat, r, values)
    if target_range == 'less than':
        arr = grad <= float(target_value)
    elif target_range == 'greater than':
        arr = grad >= float(target_value)
    #end if
    isfinite = np.isfinite(grad)
    return arr, isfinite
#end func
    
def project_onto_surface(vals, depths, min_depth=None, max_depth=None):
    """
    Function to take all depth slices of a feature between a minimum and 
    maximum depth value, and generate a single surface with all depth slices
    projected onto it.
    
    
    Parameters
    ----------
    vals: numpy array
        3D array of boolean values representing the feature at each depth 
        slice.
    
    depths: numpy array
        1D array of fepth values for each slice.
        
    min_depth: float
        Minimum depth slice to project to the surface.
        
    max_depth: float
        Maximum depth slice to project to the surface.
        
    
    Returns
    -------
    arr: numpy array
        3D numpy array (with 3rd dimension having length 1) representing the 
        surface to which the feature has been projected.
        
    
    Calls
    -----
    functions.nearest_index
    
    
    """
    if min_depth is not None:
        minind = fn.nearest_index(depths, min_depth)
    else:
        minind = 0
    #end if
    if max_depth is not None:
        maxind = fn.nearest_index(depths, max_depth)
    else:
        maxind = len(depths)
    #end if
    arr = np.expand_dims(np.any(vals[:, :, minind:maxind], axis=2), 2)
    return arr
#end func
    
def project_from_surface(vals, depths):
    """
    Function to project a feature from a surface to a set of depths. 
    
    
    Parameters
    ----------
    vals: numpy array
        3D array of boolean values representing feature which has been 
        projected to a surface. 
        
    depths: numpy array
        1D array of depth values to project the feature values to.
        
        
    Returns
    -------
    arr: numpy array
        3D numpy array representing the surfaces to which the feature has been 
        projected.
        
        
    """
    ndepths = len(depths)
    arr = np.repeat(vals, ndepths, axis=2)
    return arr
#end func

"""
def gcdist_to_feat(points, edges):
    deg_rad = np.pi/180.0
    rad_m = 6371000.0
    lon1, lat1 = points.T*deg_rad
    lon2, lat2 = edges.T*deg_rad
    colat1 = np.pi/2 - lat1
    colat2 = np.pi/2 - lat2
    min_distances = np.ones(len(points))*np.pi
    for i in range(len(edges)):
        cos1 = np.cos(colat1)
        sin1 = np.sin(colat1)
        dlon = lon1 - lon2[i]
        cos2 = np.cos(colat2[i])
        sin2 = np.sin(colat2[i])
        cos_dlon = np.cos(dlon)
        vals = cos1*cos2 + sin1*sin2*cos_dlon
        vals[vals > 1] = 1
        vals[vals < -1] = -1
        distances = np.arccos(vals)
        min_distances = np.min([min_distances, distances], axis=0)
    #end for
    return min_distances*rad_m
#end func
"""

def gcdist_to_feat(points, edges, spacing="small"):
    """
    A function to quickly calculate great circle distances from a set of points
    to the edge of the closest feature to the point.
    The cosdist_to_feat function is designed to work fast and with low memory
    usage. Since the time complexity of the algorithm used is 
    O(dx^(-2)*dy^(-2)) with dx, dy the spacing between grid points in the x and
    y directions of the geophysical model, every effort is made to reduce the
    number of calculations performed by the model.
    
    
    Parameters
    ----------
    points : numpy.ndarray
        Numpy array of shape (npoints, 2) of longitude and latitude coordinates
        of the points under consideration.
        
    edges : numpy.ndarray
        Numpy array of shape (nedges, 2) of longitude and latitude coordinates
        of the grid points which make up the edges of features in a geophysical
        model.
        
    spacing : string
        Either "small" or "large" to denote whether or not the distance between
        grid points is sufficiently large to use a less expensive formula to 
        calculate great circle distances. If "small", haversine formula is 
        used.
        
    
    Returns
    -------
    distances : numpy.ndarray
        Numpy array of shape (npoints) of distances (metres) between each point
        in points and the closest point in edges.
        
    
    Calls
    -----
    cosdist_to_feat
    
    sindist_to_feat
    
    
    """
    deg_rad = np.pi/180.0
    rad_m = 6371000.0
    lon1, lat1 = points.T*deg_rad
    lon2, lat2 = edges.T*deg_rad
    colat1 = np.pi/2 - lat1
    colat2 = np.pi/2 - lat2  
    if spacing == "large":
        cosp1 = np.cos(lon1).astype("float32")
        cost1 = np.cos(colat1).astype("float32")
        sinp1 = np.sin(lon1).astype("float32")
        sint1 = np.sin(colat1).astype("float32")
        cosp2 = np.cos(lon2).astype("float32")
        cost2 = np.cos(colat2).astype("float32")
        sinp2 = np.sin(lon2).astype("float32")
        sint2 = np.sin(colat2).astype("float32")
        in_array = -1*np.ones(len(points)).astype("float32")
        out_array = cosdist_to_feat(cosp1, cost1, sinp1, sint1, cosp2, cost2, 
                                    sinp2, sint2, in_array, in_array)
        out_array = out_array.astype("float64")
        out_array[out_array > 1] = 1
        distances = np.arccos(out_array)*rad_m
    elif spacing == "small":
        sint1 = np.sin(colat1).astype("float32")
        sint2 = np.sin(colat2).astype("float32")
        cost1_2 = np.cos(colat1/2).astype("float32")
        cost2_2 = np.cos(colat2/2).astype("float32")
        sint1_2 = np.sin(colat1/2).astype("float32")
        sint2_2 = np.sin(colat2/2).astype("float32")
        cosp1_2 = np.cos(lon1/2).astype("float32")
        cosp2_2 = np.cos(lon2/2).astype("float32")
        sinp1_2 = np.sin(lon1/2).astype("float32")
        sinp2_2 = np.sin(lon2/2).astype("float32")
        in_array = np.ones(len(points)).astype("float32")
        out_array = sindist_to_feat(sint1, sint2, cost1_2, sint1_2, cosp1_2, 
                                    sinp1_2, cost2_2, sint2_2, cosp2_2, sinp2_2, 
                                    in_array, in_array)
        out_array = out_array.astype("float64")
        distances = 2*np.arcsin(np.sqrt(out_array))*rad_m
    else:
        distances = np.zeros(len(points))
    #end if
    return distances
#end func
    
#@numba.jit(nopython=True, parallel=True)
def cosdist_to_feat(cosp1, cost1, sinp1, sint1, cosp2, cost2, sinp2, sint2, 
                    arr, d):
    for i in range(len(cosp2)):
        d = cost1*cost2[i] + sint1*sint2[i]*(cosp1*cosp2[i] + sinp1*sinp2[i])
        arr[d > arr] = d[d > arr]
    #end for
    return arr
#end func
    
#@numba.jit(nopython=True, parallel=True)
def sindist_to_feat(sint1, sint2, cost1_2, sint1_2, cosp1_2, sinp1_2, cost2_2, 
                    sint2_2, cosp2_2, sinp2_2, arr, d):
    for i in range(len(sint2)):
        d = (cost1_2*sint2_2[i] - sint1_2*cost2_2[i])**2 + \
            (sinp1_2*cosp2_2[i] - cosp1_2*sinp2_2[i])**2*sint1*sint2[i]
        arr[d < arr] = d[d < arr]
    #end for
    return arr
#end func
