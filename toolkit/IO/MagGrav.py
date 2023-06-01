import numpy as np
import pandas as pd
import toolkit.functions as fn

class MagGrav():   
    def __init__(self, model_fn, dims, cellsize, nzlst, dzlst, epsg=0, 
                 proj_str=None, bottomleft=[0, 0, 0], nanval=np.nan):
        """
        A class to facilitate the importation of gravity or magnetic susceptibility
        models, specifically in the format of the NAC data set.
        
        
        Calls
        -----
        functions.epsg_project
        
        
        """
        
        x0, y0, z0 = bottomleft
        nx, ny = dims
        nz = np.sum(nzlst)
        dx, dy = cellsize
        
        xn = x0 + nx*dx
        yn = y0 + ny*dy
        
        self.gcx = np.arange(x0, xn, dx)
        self.gcy = np.arange(y0, yn, dy)
        
        gcz = []
        for i in range(len(nzlst)):
            z = np.arange(0, nzlst[i]*dzlst[i], dzlst[i]) + z0
            z0 = z0 + nzlst[i]*dzlst[i]
            gcz.append(z)
        #end for
        self.gcz = np.hstack(gcz)
        
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
