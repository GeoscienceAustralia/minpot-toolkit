import numpy as np
import pandas as pd
import toolkit.functions as fn

class CSV_model():   
    def __init__(self, model_fn, dims, cellsize, depth_interval, epsg=0, 
                 proj_str=None, bottomleft=[0, 0, 0], nanval=np.nan):
        """
        A class to facilitate the importation of models in .csv format.
        
        
        Calls
        -----
        functions.epsg_project
        
        
        """
        
        x0, y0, z0 = bottomleft
        nx, ny, nz = dims
        dx, dy = cellsize
        dz = depth_interval
        
        xn = x0 + nx*dx
        yn = y0 + ny*dy
        zn = z0 + nz*dz
        
        self.gcx = np.arange(x0, xn, dx)
        self.gcy = np.arange(y0, yn, dy)
        self.gcz = np.arange(z0, zn, dz)
        
        df = pd.read_csv(model_fn, header=None)
        self.values = np.reshape(np.array(df), (ny, nx, nz))
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
