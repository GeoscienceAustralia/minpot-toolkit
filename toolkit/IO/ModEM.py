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

import os
import toolkit.functions as fn
from mtpy.modeling.modem import Data, Model
import numpy as np
from shapely.geometry import MultiPoint, Point

class ModEM():
    def __init__(self,model_fn,data_fn=None,stationxy_fn=None,proj_str=None,
                 epsg=None,rpad=0,centre=None,stationfile_kwargs={},
                 stationll_fn=None):
        """
        A class to facilitate the importation of ModEM geophysical models.
 
        
        Calls
        -----
        mtpy
        
        functions.epsg_project
        
        
        """
        self.wd = os.path.dirname(model_fn)
        
        self.mObj = Model()
        self.mObj.read_model_file(model_fn=model_fn)
        
        if data_fn is not None:
            self.dObj = Data(model_epsg=epsg)
            self.dObj.read_data_file(data_fn=data_fn)
            self.slon = self.dObj.station_locations.lon
            self.slat = self.dObj.station_locations.lat
            
            
            if centre is None:
                if epsg == 0:
                    centre = fn.epsg_project(self.dObj.center_point['lon'],
                                             self.dObj.center_point['lat'],
                                             4326, 0, proj_str=proj_str)
                else:
                    centre = [self.dObj.center_point['east'],
                                  self.dObj.center_point['north']]
        else:
            self.dObj = None
            if stationll_fn is not None:
                self.slon,self.slat = np.genfromtxt(stationll_fn,unpack=True,
                                                 **stationfile_kwargs)
            if stationxy_fn is not None:
                sx,sy = np.genfromtxt(stationxy_fn,unpack=True,
                                                 **stationfile_kwargs)
                self.slon,self.slat = fn.epsg_project(sx+centre[0],
                                                      sy+centre[1],
                                                      0,4326,proj_str=proj_str)
            
        
        # wrap station longitude
        
        if np.any(self.slon > 180):
            gt180 = np.where(self.slon > 180)[0]
            
            self.slon[gt180] = self.slon[gt180] % 360 - 360
        
        if epsg is None:
            epsg=0
        self.epsg = epsg
        
        self.gx = self.mObj.grid_east + centre[0]
        self.gy = self.mObj.grid_north + centre[1]
        self.gz = self.mObj.grid_z
        
        i0,i1 = rpad, self.gx.size - rpad
        j0,j1 = rpad, self.gy.size - rpad
        
        self.gx = self.gx[i0:i1]
        self.gy = self.gy[j0:j1]
        
        self.gx_grid, self.gy_grid = np.meshgrid(self.gx, self.gy)
        
        self.gcx, self.gcy, self.gcz = [np.array([np.mean(arr[i:i+2]) for i in\
                                        range(len(arr)-1)]) for arr in \
                                        [self.gx,self.gy,self.gz]]
        self.gcx_grid, self.gcy_grid = np.meshgrid(self.gcx, self.gcy)

        self.glon, self.glat = fn.epsg_project(self.gx_grid, self.gy_grid,
                                               self.epsg, 4326,
                                               proj_str=proj_str)
        self.gclon, self.gclat = fn.epsg_project(self.gcx_grid, self.gcy_grid,
                                                 self.epsg, 4326,
                                                 proj_str=proj_str)
        self.sx, self.sy = fn.epsg_project(self.slon,self.slat,4326,self.epsg,
                                           proj_str=proj_str)



        
        self.values = self.mObj.res_model[j0:j1-1,i0:i1-1]
        
    def mask_station_distance(self,distance_degrees):
        
        # station locations
        station_multipoint = MultiPoint(np.vstack([self.slon,self.slat]).T)
        
        # mask within distance_degrees of a station
        self.resmask = np.zeros_like(self.values[:,:,0],dtype = bool)
        for ii in range(len(self.gclon)):
            for jj in range(len(self.gclat[ii])):
                if station_multipoint.distance(Point(self.gclon[ii,jj],
                                                     self.gclat[ii,jj])) > distance_degrees:
                    self.resmask[ii,jj] = True
                    
        self.values[self.resmask] = np.nan
        self.bounds = [np.min(self.slon) - distance_degrees, 
                       np.max(self.slon) + distance_degrees, 
                       np.min(self.slat) - distance_degrees, 
                       np.max(self.slat) + distance_degrees]
        
        
    def write_station_xy(self):
        np.savetxt(os.path.join(self.wd, 'Station_XY.txt'),
                   np.vstack([self.dObj.station_locations.rel_east,
                              self.dObj.station_locations.rel_north]).T,
                   fmt='%.2f')
        
   
    def write_station_ll(self):
        np.savetxt(os.path.join(self.wd, 'Station_LL.txt'),
                   np.vstack([self.slon,
                              self.slat]).T,
                   fmt='%.6f')
#end class
