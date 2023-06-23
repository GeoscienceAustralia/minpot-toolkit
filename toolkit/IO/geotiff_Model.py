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
import toolkit.functions as fn
import rasterio
from affine import Affine
from pyproj import Proj, transform

class geotiff_model():
    def __init__(self, model_fn, interval=1, force_epsg=False):
        """
        A class to facilitate the importation of models in .tif format.


        Calls
        -----
        functions.epsg_project


        """
        rh = rasterio.open(model_fn)
        epsg = rh.crs.to_epsg()
        T0 = rh.transform
        p1 = Proj(rh.crs)
        td = rh.read()[0, :, :]

        cols, rows = np.meshgrid(np.arange(td.shape[1]), np.arange(td.shape[0]))
        T1 = T0 * Affine.translation(0.5, 0.5)

        x, y = T1*(cols.flatten(), rows.flatten())
        xgrid = np.reshape(x, cols.shape)
        ygrid = np.reshape(y, rows.shape)

        # drop samples
        xgrid = xgrid[::interval, ::interval]
        ygrid = ygrid[::interval, ::interval]
        td = td[::interval, ::interval]

        self.gcx = xgrid[0,:]
        self.gcy = ygrid[:,0]
        self.gcz = np.array([0.0])
        self.values = np.expand_dims(td, 2)

        if self.gcy[1] - self.gcy[0] < 0:
            self.gcy = self.gcy[::-1]
            self.values = self.values[::-1, :]
            ygrid = ygrid[::-1, :]
        #end if

        if force_epsg:
            self.epsg = force_epsg
        else:
            self.epsg = epsg
        #end if
        self.gclon, self.gclat = fn.epsg_project(xgrid, ygrid, self.epsg, 4326)

        for nodataval in rh.nodatavals:
            self.values[self.values == nodataval] = np.nan
        # end for

        minLon = np.min(self.gclon)
        minLat = np.min(self.gclat)
        maxLon = np.max(self.gclon)
        maxLat = np.max(self.gclat)
        self.bounds = [minLon, maxLon, minLat, maxLat]
    #end func
#end class
