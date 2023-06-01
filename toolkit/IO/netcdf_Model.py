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
