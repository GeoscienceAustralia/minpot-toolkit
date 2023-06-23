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
Prospectivity module
-----------------
This module is used to produce an area reduction metric for the search 
space when prospecting for minerals, based on the results of the statistics 
module.
    
"""

import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.stats import norm
import toolkit.plotting as pt

class Prospectivity():
    def __init__(self, stat_obj_dict, feature_dict):
        """
        A class to facilitate analysing the prospectivity of an area based on
        the results of statistical tests performed.
        
        
        Parameters
        ----------
        stat_obj_dict: dictionary
            Dictionary containing an object for each statistical test which was
            performed.
            
        
        Calls
        -----
        prospectivity_map
        
        
        """
        keys = stat_obj_dict.keys()
        self.dct = {}
        for key in keys:
            if stat_obj_dict[key].prospectivity_map:
                feature_key = stat_obj_dict[key].feature_key
                self.dct[key] = prospectivity_map(stat_obj_dict[key],
                                                  feature_dict[feature_key])
                self.dct[key].Process()
            #end if
        #end for
    #end func
    
    def Plots(self, bgimage="", wd_images='.'):
        """
        Generates a prospectivity map based on the results of each statistical
        test that has been performed.
        
        
        Parameters
        ----------
        bgimage: string or matplotlib.imread object
            Background image to overlay the generated map onto.
            
        wd_images: string
            Directory in which to save figure.


        """
        for key in self.dct.keys():
            self.dct[key].generate_maps(bgimage=bgimage, wd_images=wd_images)
        #end for
    #end func
#end class
    
class prospectivity_map():
    def __init__(self, stat_obj, feature):
        """
        A class to produce a prospectivity map based on input criteria.
        
        
        Parameters
        ----------
        stat_obj: stat_corr (various) object
            Object containing the results of a statistical test which has been 
            performed.
            
        
        """
        self.stat_obj = stat_obj
        self.lon = feature.gclon
        self.lat = feature.gclat
        self.distances = feature.distances
    #end func
    
    def compute_and_interpolate_quantiles(self):
        """
        Finds the quantiles (in steps of 1%) for distances from deposit 
        locations to the feature under investigation.
        
        
        Calls
        -----
        compute_quantiles
        
        interpolate_quantiles
        
        
        """
        if not hasattr(self, 'distances'):
            self.compute_distances_to_feature()
        #end if
        if not hasattr(self, 'dist_arr'):
            self.dist_arr = self.stat_obj.distance_array['distance'][0]
        #end if
        dq = 0.01
        quantiles = np.arange(0, 1+dq, dq)
        self.quantiles = compute_quantiles(self.dist_arr, quantiles)
        self.fun, self.funinv = interpolate_quantiles(self.quantiles)
    #end func
    
    def find_area_enclosed_by_quantile(self):
        """
        Function to calculate the area which is enclosed by a distance 
        quantile.
        
        
        """
        if not hasattr(self, 'quantiles'):
            self.compute_and_interpolate_quantiles()
        #end if
        prop = list()
        dims = self.distances.shape
        n_grid_points = dims[0]*dims[1]*dims[2]
        for row in self.quantiles:
            d = row[1]
            arr = (self.distances <= d)
            num = np.sum(arr)
            prop.append(num/n_grid_points)
        #end for
        self.prop = np.array(prop)
    #end func
    
    def compute_and_interpolate_pdf(self):
        """
        Computes and interpolates the probability density function of the 
        distances to the feature from mineral deposit locations.
        
        
        Calls
        -----
        compute_and_interpolate_pdf
        
        
        """
        if not hasattr(self, 'dist_arr'):
            self.dist_arr = self.stat_obj.distance_array['distance'][0]
        #end if
        self.pdf = compute_and_interpolate_pdf(self.dist_arr, smooth=True)
    #end if
    
    def generate_maps(self, bgimage="", wd_images='.'):
        """
        Creates a figure showing the quantile each grid point in an array falls
        within.
        
        
        Parameters
        ----------
        bgimage: string or matplotlib.imread object
            Background image to overlay the map onto.
            
        wd_images: string
            Directory in which to save figure.
            
        
        Calls
        -----
        plotting.plot_quantile_map
        
        plotting.plot_pdf_map
        
        
        """
        if not hasattr(self, 'fun'):
            self.compute_and_interpolate_quantiles()
        #end if        
        wd_images = os.path.join(wd_images, 'Stats')
        self.quantile_values = 1 - self.fun(self.distances)
        pt.plot_quantile_map(self, bgimage=bgimage, wd_images=wd_images)
        self.pdf_values = self.pdf(self.distances)
        pt.plot_pdf_map(self, bgimage=bgimage, wd_images=wd_images)
    #end func
    
    def Process(self):
        """
        Computes distances to a feature from each grid point in an array, and 
        finds the quantiles (in steps of 1%) for distances from deposit 
        locations to the feature under investigation.
        
        
        """
        self.compute_and_interpolate_quantiles()
        self.compute_and_interpolate_pdf()
    #end func
#end class
    
def compute_and_interpolate_pdf(dist_arr, smooth=False):
    """
    Computes a probability density function by computing a cumulative 
    distribution function and differentiating.
    
    
    Parameters
    ----------
    dist_arr : numpy array
        Array of distance values from deposits to a feature.
        
    
    Returns
    -------
    pdf : scipy.interpolate.interp1d object
        Probability density function
    
    
    """
    if not smooth:
        x = np.sort(dist_arr)
        np.save('DIST_ARRAY.npy',x)
        y = np.linspace(0, 1, len(dist_arr))
        dydx1 = np.array([(y[1] - y[0])/(x[1] - x[0])])
        dydxn = np.array([(y[-1] - y[-2])/(x[-1] - x[-2])])
        dydx = (y[2:] - y[:-2])/(x[2:] - x[:-2])
        dydx = np.hstack([dydx1, dydx, dydxn])
        pdf = interp1d(x, dydx, fill_value=0, bounds_error=False, kind='linear')
    else:
        mean = dist_arr.mean()
        std  = dist_arr.std()
        min  = dist_arr.min()
        max  = dist_arr.max()
        rng  = max-min
        X    = np.linspace(min - 0.25*rng, max + 0.25*rng, 1001)
        rv   = norm(mean, std)
        Y    = rv.pdf(X)
        pdf  = interp1d(X, Y, fill_value=0, bounds_error=False, kind='linear')
    #end if
    return pdf
#end func
    
def compute_quantiles(arr, quantiles):
    """
    Computes quantiles for an array.
    
    
    Parameters
    ----------
    arr: numpy array
        Array to compute quantiles for.
        
    quantles: numpy array
        Array containing values for each quantile we require.
        
    
    Returns
    -------
    quantiles: numpy array
        Quantiles for the array.
        
    
    """
    dvals = list()
    for q in quantiles:
        d = np.quantile(arr, q)
        dvals.append(d)
    #end for
    dvals = np.array(dvals)
    return np.vstack([quantiles, dvals]).T
#end func
    
def interpolate_quantiles(quantiles):
    """
    Takes an array of x and y coordinates, and interpolates z=1-y values.
    
    
    Parameters
    ----------
    quantiles: numpy array
        2D containing (x, y) values for x (distance) and y (quantile which
        x falls within for the deposit dataset distance array).
        
    
    Returns
    -------
    fun: scipy.interpolate.interp1d object
        A linear interpolant for the quantiles.
        
        
    """
    y, x = quantiles.T
    fun = interp1d(x, y, fill_value=(0, 1), bounds_error=False, kind='linear')
    funinv = interp1d(y, x, fill_value=(0, 1), bounds_error=False, 
                      kind='linear')
    return fun, funinv
#end func
