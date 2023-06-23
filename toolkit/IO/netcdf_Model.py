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

import numpy as np
from netCDF4 import Dataset
import toolkit.functions as fn

class netcdf_model():   
    def __init__(self, model_fn, x_var_name, y_var_name, values_var_name, 
                 z_var_name=None, epsg=0, proj_str=None, nanval=np.nan,
                 swap_xy=False, interval=1):
        """
        A class to facilitate the importation of models in .nc format.
        
        
        Calls
        -----
        functions.epsg_project
        
        
        """
        
        ds = Dataset(model_fn, mode='r')
        
        self.gcx = np.array(ds.variables[x_var_name][:])
        self.gcy = np.array(ds.variables[y_var_name][:])
        self.values = np.array(ds.variables[values_var_name][:])
        self.gcx = self.gcx[::interval]
        self.gcy = self.gcy[::interval]
        if z_var_name is not None:
            self.gcz = np.array(ds.variables[z_var_name][:])
            self.values = self.values[::interval, ::interval]
        else:
            self.gcz = np.array([0.0])
            self.values = np.expand_dims(self.values, 2)
            self.values = self.values[::interval, ::interval, :]
        #end if
        
        if swap_xy:
            self.values = self.values.transpose(1, 0, 2)
        #end if
        
        if self.gcy[1] - self.gcy[0] < 0:
            self.gcy = self.gcy[::-1]
            self.values = self.values[::-1, :]
        #end if
        
        self.values[self.values == nanval] = np.nan
        
        xgrid, ygrid = np.meshgrid(self.gcx, self.gcy)
        self.gclon, self.gclat = fn.epsg_project(xgrid, ygrid, epsg, 4326, 
                                                 proj_str)
        
        minLon = np.min(self.gclon)
        minLat = np.min(self.gclat)
        maxLon = np.max(self.gclon)
        maxLat = np.max(self.gclat)
        self.bounds = [minLon, maxLon, minLat, maxLat]
    #end func
#end class
