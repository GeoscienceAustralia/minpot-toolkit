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
Model input module
-----------------
This module is used to import geophysical models from various input formats.
Supported formats are comma separated values (.csv), ModEM, and netcdf (.nc).
    
"""

import os, time
import cartopy.io.shapereader as shpreader
import numpy as np
import toolkit.functions as fn

from geopandas import read_file
from geopandas.geoseries import GeoSeries
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.ops import unary_union

class ModelsManager():
    def __init__(self):
        """
        Class to faciliate importation of geophysical models.
        
        
        """
        self.cellsize = 0.1
        self.bounds = np.array([110.0, 155.0, -45.0, -10.0])
        self.outside_polygon_file = ""
        self.inside_polygon_file = ""
    #end func
    
    def compute_domain(self):
        """
        Takes polygons of region of coverage for each model, merges them, and 
        computes the intersection of this domain with land.
        
        
        Calls
        -----
        get_outside_polygon
        
        
        """
        t0 = time.time()
        print('Generating land polygon')
        land_shp_fname = \
            shpreader.natural_earth(resolution='50m', category='physical', 
                                    name='land')
        land_geom = \
            unary_union(list(shpreader.Reader(land_shp_fname).geometries()))
        t1 = time.time()
        print('Finished generating land polygon, time =', t1-t0)
        print('Merging model domains')
        domain = unary_union(self.domain)
        xmin, ymin, xmax, ymax = domain.bounds
        self.bbox = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), 
                             (xmin, ymax)])
        t2 = time.time()
        print('Finished merging model domains, time =', t2-t1)
        print('Computing intersection of domain with specified shapefiles')
        land = self.bbox.intersection(land_geom)
        region = domain.intersection(land)        
        self.domain = {'region': region}
        if self.outside_polygon_file:
            outside_polygon = get_polygon(self.outside_polygon_file)
            allowed_region = region.difference(outside_polygon)
        else:
            allowed_region = region
        #end if
        if self.inside_polygon_file:
            inside_polygon = get_polygon(self.inside_polygon_file)
            allowed_region = allowed_region.intersection(inside_polygon)
        #end if
        self.domain['allowed_region'] = allowed_region
        t3 = time.time()
        print('Completed computing region, time =', t3-t2)
    #end func
    
    def decompose_domain_into_grid(self):
        """
        Takes the region of coverage of the model, and fishnets the polygon to
        improve the speed of the filtering of random points by location.
        
        
        """
        print('Decomposing region into grid')
        t3 = time.time()
        edge_squares = [self.bbox]
        gs_within = GeoSeries()
        allowed_region = self.domain['allowed_region']
        not_allowed_region = self.bbox.difference(allowed_region)
        i = 0
        t4 = time.time()
        for d in [1, 0.5, 0.25]:
            i = i + 1
            dx = dy = d
            squares = list()
            for edge_square in edge_squares:
                xmin, ymin, xmax, ymax = edge_square.bounds
                xvals = np.arange(xmin, xmax+dx, dx)
                yvals = np.arange(ymin, ymax+dy, dy)
                squares.extend([Polygon([(x,y), (x+dx,y), (x+dx,y+dy), 
                                         (x,y+dy)]) \
                                for x in xvals for y in yvals])
            #end for
            gs = GeoSeries(squares)
            ind1 = gs.intersects(allowed_region)
            ind2 = gs.intersects(not_allowed_region)
            edge_filt = np.logical_and(ind1, ind2)
            within_filt = np.logical_and(ind1, ~ind2)
            gs_within = gs_within.append(gs[within_filt])
            edge_squares = list(gs[edge_filt])
            print('Iteration', i, 'complete of', 3, 'time =', time.time()-t4)
            t4 = time.time()
        #end for
        print('Completed generating grid, time =', t4-t3)
        print('Computing union of grid squares')
        gs = gs_within.append(GeoSeries(edge_squares))
        grid = unary_union(gs)
        t5 = time.time()
        print('Completed computing union of grid squares, time =', t5-t4)
        self.domain['grid'] = grid
    #end func
    
    def ParseXMLNode(self, rootNode, commandList):
        depths_type = rootNode.attrib["depths_type"]
        if depths_type == "values":
            self.depth_list = fn.changetype(rootNode.attrib["depths"], 'array', 
                                            dtype=float)
        elif depths_type == "range":
            ranges = fn.changetype(rootNode.attrib["depths"], list, dtype=str, 
                                   sep=';')
            depths = list()
            for rng in ranges:
                minval, maxval, interval = \
                    fn.changetype(rng, 'array', dtype=float)
                depths.append(np.arange(minval, maxval + interval, interval))
            #end for
            self.depth_list = np.hstack(depths)
        #end if
        models = rootNode.findall("model")
        for row in commandList:
            if hasattr(self, row[0]):
                tp = type(getattr(self, row[0]))
                val = fn.changetype(row[1], tp)
                setattr(self, row[0], val)
            #end if
        #end for
        model_details = list()
        for model in models:
            name = model.attrib["name"]
            model_type = model.attrib["type"]
            model_details.append([str(name + '.' + model_type), model])
        #end for
        self.model_details_dict = {row[0]: row[1] for row in model_details}
    #end func
    
    def Read_Models(self):
        """
        Reads in all input models, creating a different geophysical model for 
        each model type.
        
        Calls
        -----
        GeophysicalModel
        
        
        """
        self.model_dict = {}
        self.domain = []
        self.model_types = [key.split('.')[-1] for key in 
                            self.model_details_dict]
        for model_type in self.model_types:
            keys = [key for key in self.model_details_dict if \
                    key.split('.')[-1] == model_type]
            models = [self.model_details_dict[key] for key in keys]
            GM = GeophysicalModel(models, model_type, self.depth_list,
                                  cellsize = self.cellsize,
                                  bounds = self.bounds)
            GM.merge_models()
            self.domain.append(GM.domain)
            self.model_dict[model_type] = \
                {'longitude': GM.longitude, 'latitude': GM.latitude,
                 'depth': GM.depth_list, 'values': GM.values,
                 'propertyname': GM.propertyname}
        #end for
    #end func
    
    def Process(self):
        self.Read_Models()
        self.compute_domain()
        self.decompose_domain_into_grid()
    #end func
#end class

class GeophysicalModel():
    def __init__(self, models, model_type, depth_list, **kwargs):
        """
        Class to facilitate importation of a geophysical model.
        
        
        Parameters
        ----------
        models : dictionary
            Dictionary of model details describing various parameters, paths,
            file names, etc.
            
        model_type : string
            Keyword for the type of model to be generated. Used to extract the
            correct information from the models dictionary.
            

        """
        
        self.model_type = model_type
        self.models = models
        self.propertyname = None
        self.cellsize = 0.1
        self.depth_list = depth_list
        self.bounds = [110, 155, -45, -10]
        self.epsg = False
        
        for key in kwargs.keys():
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
            #end if
        #end for

    def read_models(self):
        """
        Determines type of model, and reads inputs into a dictionary.


        Calls
        -----
        read_modem
        
        read_maggrav
        
        read_csv_model
        
        read_netcdf_model
        
        
        """
        self.model_dict = {}
        if self.model_type == 'ModEM':
            self.propertyname = 'resistivity'
        elif self.model_type == 'Grav':
            self.propertyname = 'density'
        elif self.model_type == 'Mag':
            self.propertyname = 'susceptibility'
        else:
            self.propertyname = 'unknown'
        #end if
        for model in self.models:
            name = model.attrib['name']
            if model.attrib['model_format'] == 'ModEM':
                self.model_dict[name] = read_modem(model)
            elif model.attrib['model_format'] == 'MagGrav':
                self.model_dict[name] = read_maggrav(model)
            elif model.attrib['model_format'] == 'csv':
                self.model_dict[name] = read_csv_model(model)
            elif model.attrib['model_format'] == 'netcdf':
                self.model_dict[name] = read_netcdf_model(model)
            elif model.attrib['model_format'] == 'geotiff':
                print(name)
                self.model_dict[name] = read_geotiff_model(model)
            # Correct for longitude wrap-around
            if np.any(self.model_dict[name].gclon > 180):
                # Extract subset for longitude wrapping
                wrapBool = self.model_dict[name].gclon > 180
                boolShape    = (wrapBool.shape[0], wrapBool[0].sum())
                invBoolShape = (wrapBool.shape[0], (wrapBool==False)[0].sum())
                wrapLon  = self.model_dict[name].gclon[wrapBool] - 360
                wrapLon = wrapLon.reshape(boolShape)
                # Re-arange model arrays
                gclon  = np.hstack((wrapLon, 
                                    self.model_dict[name].gclon[wrapBool==False].reshape(invBoolShape)))
                gclat  = np.hstack((self.model_dict[name].gclat[wrapBool].reshape(boolShape), 
                                    self.model_dict[name].gclat[wrapBool==False].reshape(invBoolShape)))
                values = np.hstack((self.model_dict[name].values[wrapBool].reshape(boolShape), 
                                    self.model_dict[name].values[wrapBool==False].reshape(invBoolShape)))
                self.model_dict[name].gclon  = gclon
                self.model_dict[name].gclat  = gclat
                self.model_dict[name].values  = np.expand_dims(values, 2)
                gcxBool = (self.model_dict[name].gcx >= 180) & (self.model_dict[name].gcx < 360)
                self.model_dict[name].gcx = np.hstack((self.model_dict[name].gcx[gcxBool] - 360,
                                                       self.model_dict[name].gcx[self.model_dict[name].gcx <= 180]))
                # Update model bounds
                minLon = np.min(self.model_dict[name].gclon)
                minLat = np.min(self.model_dict[name].gclat)
                maxLon = np.max(self.model_dict[name].gclon)
                maxLat = np.max(self.model_dict[name].gclat)
                self.model_dict[name].bounds = [minLon, maxLon, minLat, maxLat]
            #end if
        #end for
    #end func
    
    def compute_coverage_region(self):
        """
        Creates a multipolygon of a triangulation of the region of coverage
        for each model.
        
        
        Calls
        -----
        compute_coverage_polygon
        
        
        """

    #end func
    
    def interpolate_to_grid(self):
        """
        Takes separate models and placed them all on the same grid.
        
        
        Calls
        -----
        interpolate_to_grid
        
        
        """
        print('Interpolating to grid')
        log = (self.model_type == 'ModEM')
        self.longitude, self.latitude, self.values = \
            interpolate_to_grid(self.cellsize, self.depth_list, 
                                self.model_dict, self.model_type, self.bounds, 
                                log=log)
        self.gclon, self.gclat = np.meshgrid(self.longitude, self.latitude)
        print('Finished interpolating to grid')
    #end func
    
    def merge_models(self):
        """
        Reads in all data, computes a multipolygon of the region of coverage,
        and places all models on the same grid.
        
        
        """
        self.read_models()
        self.interpolate_to_grid()
        print('Computing coverage region')
        self.domain = compute_coverage_polygon(self)
        print('Completed computing coverage region')
    #end func
#end class
    
def compute_coverage_polygon(model):
    """
    Function to generate a polygon of the region of coverage for a model.
    Note - if the longitude values cross over the 180 degree line, then the 
    furthest eastern cells in the eastern hemisphere are not included in the 
    region of coverage.
    
    
    Parameters
    ----------
    model_dict : dictionary
        Dictionary of model objects to use. Must contain longitude and latitude
        values for grid points in a mesh grid, and function values in a rank-3
        tensor (3rd dimension may be 1).
        
    
    Returns
    -------
    domain : shapely multipolygon
        Multipolygon of the region of coverage for the model.
        
        
    Calls
    -----
    functions.gradient
    
    
    """
    
    bounds = model.bounds
    gclon = model.gclon
    gclat = model.gclat
    values = model.values
    arr = np.any(~np.isnan(values), axis=2)
    grad = (fn.gradient(gclon, gclat, np.array([0]),
                        1.0*np.expand_dims(arr, 2))[:,:,0] != 0)*1
    iindices = np.arange(len(grad[:, 0]))
    jindices = np.arange(len(grad[0, :]))
    grad1 = np.zeros_like(grad)
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            grad1 = grad1 + \
                grad.take(iindices+i, axis=0, 
                          mode='wrap').take(jindices+j, axis=1, 
                                            mode='wrap')
        #end for
    #end for
    on_boundary = np.zeros_like(grad1).astype(bool)
    on_boundary[:, 0] = True
    on_boundary[:, -1] = True
    on_boundary[0, :] = True
    on_boundary[-1, :] = True
    near_edge = (grad1 > 0)
    not_inside_edges = np.logical_or((grad > 0), ~arr)
    vals = np.logical_or(np.logical_and(near_edge, not_inside_edges),
                         on_boundary)
    
    xfilt = gclon[vals]
    yfilt = gclat[vals]
    zfilt = arr[vals]
    points = np.transpose(np.vstack([xfilt, yfilt]))
    tri = Delaunay(points)
    indices = [np.all(zfilt[simp]) for simp in tri.simplices]
    coord_groups = [tri.points[simp] for simp in tri.simplices[indices]]
    polygons = [Polygon(coords) for coords in coord_groups]
    
    region = unary_union(polygons)
    minlon, maxlon, minlat, maxlat = bounds
    corners = [(minlon, minlat), (minlon, maxlat), (maxlon, maxlat), 
           (maxlon, minlat)]
    polygon = Polygon(corners)
    
    return region.intersection(polygon)
#end func
    
def get_polygon(filename):
    shapefile = read_file(filename)
    return unary_union(shapefile['geometry'])
#end func
    
def interpolate_to_grid(cellsize, depth_list, model_dict, model_type, bounds, 
                        log=True):
    """
    Interpolate values from a geophysical model to a grid.
    

    Parameters
    ----------
    cellsize : float
        Width of cells for grid.
        
    depth_list : list of float
        List of depth values data is provided at.
        
    model_dict : dictionary 
        Dictionary of geophysical models.
        
    model_type : string
        String identifying the type of model being used.
        
    log : boolean
        Describes whether to interpolate function values or their logarithm.
        
        
    Returns
    -------
    longitude : array
        Numpy linspace of longitude values in grid.
    
    latitude : array
        Numpy linspace for latitude values in grid.
    
    values : array
        Rank-3 tensor of interpolated values for points on the grid.
        
        
    Calls
    -----
    functions.epsg_project
    
        
    """
    
    # make a master grid to store data in
    #minLon, minLat, maxLon, maxLat = -180, -90, 180, 90
    globMinLon, globMaxLon, globMinLat, globMaxLat = bounds
    if globMaxLon < globMinLon:
        globMaxLon = globMaxLon + 360
    #end if
    
    longitude = np.linspace(globMinLon, globMaxLon, 
                            int((globMaxLon-globMinLon)/cellsize+1))
    latitude = np.linspace(globMinLat, globMaxLat, 
                           int((globMaxLat-globMinLat)/cellsize+1))
    values = np.zeros((longitude.size, latitude.size, depth_list.size))
    values[:] = np.nan
    
    # dictionary to store min/max lat/lon in each model
    for key in model_dict.keys():
        modObj = model_dict[key]
        if modObj.epsg != 0:
            epsg = modObj.epsg
            proj_str = None
        else:
            epsg = 0
            proj_str = modObj.proj_str
        #end if
    
        # get min/max longitude and latitude
        minLon, maxLon, minLat, maxLat = modObj.bounds
    
        # force onto even lat/lon intervals
        minLon, minLat = [np.floor(val) for val in [minLon, minLat]]
        maxLon, maxLat = [np.ceil(val) for val in [maxLon, maxLat]]
        if maxLon <= minLon:
            maxLon = maxLon + 360
        #end if
        
        minLon = max(minLon, globMinLon)
        minLat = max(minLat, globMinLat)
        maxLon = min(maxLon, globMaxLon)
        maxLat = min(maxLat, globMaxLat)

        # define longitude and latitude to interpolate to
        loni = np.linspace(minLon, maxLon, int((maxLon-minLon)/cellsize + 1))
        lati = np.linspace(minLat, maxLat, int((maxLat-minLat)/cellsize + 1))
        loni, lati, zi = np.meshgrid(loni, lati, depth_list)   
        
        #Assumption is made that x, y coordinates do not wrap around in
        #original model
        xi, yi = fn.epsg_project(loni, lati, 4326, epsg, proj_str=proj_str)

        if log:
            vals = np.log10(modObj.values)
        else:
            vals = modObj.values
        #end if
        
        if modObj.gcz.shape[0] > 1:
            func = RegularGridInterpolator((modObj.gcy, modObj.gcx, 
                                            modObj.gcz), vals, 
                                           bounds_error=False)
            
            xyzi = np.vstack([yi.flatten(), xi.flatten(), zi.flatten()]).T
            ny, nx, nz = xi.shape
            
            if log:
                newvals = 10**func((xyzi)).reshape(ny, nx, nz).transpose(1,0,2)
            else:
                newvals = func((xyzi)).reshape(ny, nx, nz).transpose(1,0,2)
            #end if
        else:
            func = RegularGridInterpolator((modObj.gcy, modObj.gcx), 
                                           vals[:,:,0], bounds_error=False)
            
            xyi = np.vstack([yi.flatten(), xi.flatten()]).T
            ny, nx, nz = xi.shape
            
            if log:
                newvals = 10**func((xyi)).reshape(ny, nx, nz).transpose(1,0,2)
            else:
                newvals = func((xyi)).reshape(ny, nx, nz).transpose(1,0,2)
            #end if
        #end if            

        # mesh into master grid
        
        i0 = fn.nearest_index(loni[0, 0, 0], longitude)
        i1 = i0 + loni.shape[1]
        
        j0 = fn.nearest_index(lati[0, 0, 0], latitude)
        j1 = j0 + lati.shape[0]    
        
        if i1 <= len(longitude):
            values[i0:i1, j0:j1][np.isfinite(newvals)] = \
                newvals[np.isfinite(newvals)]
        else:
            newvals_east = newvals[:len(longitude)-i0, :, :]
            newvals_west = newvals[len(longitude)-i1:, :, :]
            values[i0:, j0:j1][np.isfinite(newvals_east)] = \
                newvals_east[np.isfinite(newvals_east)]
            values[:i1-len(longitude), j0:j1][np.isfinite(newvals_west)] = \
                newvals_west[np.isfinite(newvals_west)]
        #end if
    #end for
    
    while np.any(longitude[longitude > 180]):
        longitude[longitude > 180] = longitude[longitude > 180] - 360
    
    values = values.transpose(1,0,2)

    return longitude, latitude, values
#end func
    
def read_csv_model(model):
    """
    Read CSV models and return CSV_model object.

    
    Parameters
    ----------
    model : xml node
        Node containing model details.
        
        
    Returns
    -------
    modObj : CSV_model object.
        
        
    Calls
    -----
    IO.CSV_Model.CSV_model
    
    """
    from toolkit.IO.CSV_Model import CSV_model
    
    wd = model.attrib["path"]
    model_fn = model.attrib["model_file_name"]
    dims = fn.changetype(model.attrib["dims"], "array", dtype=int)
    cellsize = fn.changetype(model.attrib["cellsize"], "array", dtype=float)
    depth_interval = fn.changetype(model.attrib["depth_interval", float])
    bottomleft = fn.changetype(model.attrib["bottomleft"], "array", 
                               dtype=float)
    nanval = fn.changetype(model.attrib["nanval"], float)
    if "epsg" in model.attrib:
        proj_str = None
        epsg = fn.changetype(model.attrib["epsg"], int)
    else:
        proj_str = model.attrib["proj_str"]
        epsg = 0
    #end if
    
    modObj = CSV_model(os.path.join(wd, model_fn), dims, cellsize, 
                       depth_interval, epsg=epsg, proj_str=proj_str, 
                       bottomleft=bottomleft, nanval=nanval)
    
    modObj.epsg = epsg
    modObj.proj_str = proj_str
    return modObj
#end func
    
def read_maggrav(model):
    """
    Read MagGrav models and return MagGrav object.

    
    Parameters
    ----------
    model : xml node
        Node containing model details.
        
        
    Returns
    -------
    modObj : MagGrav object.
        
        
    Calls
    -----
    IO.MagGrav.MagGrav
    
    """
    from toolkit.IO.MagGrav import MagGrav
    
    wd = model.attrib["path"]
    model_fn = model.attrib["model_file_name"]
    dims = fn.changetype(model.attrib["dims"], "array", dtype=int)
    cellsize = fn.changetype(model.attrib["cellsize"], "array", dtype=float)
    bottomleft = fn.changetype(model.attrib["bottomleft"], "array", 
                               dtype=float)
    nanval = fn.changetype(model.attrib["nanval"], float)
    nzlst = fn.changetype(model.attrib["nzlst"], "array", dtype=int)
    dzlst = fn.changetype(model.attrib["dzlst"], "array", dtype=float)
    
    if "epsg" in model.attrib:
        proj_str = None
        epsg = fn.changetype(model.attrib["epsg"], int)
    else:
        proj_str = model.attrib["proj_str"]
        epsg = 0
    #end if
    
    modObj = MagGrav(os.path.join(wd, model_fn), dims, cellsize, nzlst, dzlst, 
                     epsg=epsg, proj_str=proj_str, bottomleft=bottomleft, 
                     nanval=nanval)
    
    modObj.epsg = epsg
    modObj.proj_str = proj_str
    return modObj
#end func

def read_modem(model):
    """
    Read ModEM models and return ModEM object.

    
    Parameters
    ----------
    model : xml node
        Node containing model details.
        
        
    Returns
    -------
    modObj : ModEM object.
        
        
    Calls
    -----
    IO.ModEM.ModEM
    
    """
    from toolkit.IO.ModEM import ModEM
    
    wd = model.attrib["path"]
    model_fn = model.attrib["model_file_name"]
    stationxy_fn = model.attrib["stationxy_fn"]
    centre = fn.changetype(model.attrib["centre"], "array", dtype=float)
    station_buffer = fn.changetype(model.attrib["station_buffer"], float)
    
    if "epsg" in model.attrib:
        proj_str = None
        epsg = fn.changetype(model.attrib["epsg"], int)
    else:
        proj_str = model.attrib["proj_str"]
        epsg = 0
    #end if
    
    modObj = ModEM(os.path.join(wd, model_fn),
                   stationxy_fn=os.path.join(wd, stationxy_fn),
                   proj_str=proj_str, epsg=epsg, rpad=1, centre=centre)
    
    modObj.mask_station_distance(station_buffer)
    modObj.proj_str = proj_str
    return modObj
#end func
    
def read_netcdf_model(model):
    """
    Read netcdf models and return netcdf_model object.

    
    Parameters
    ----------
    model : xml node
        Node containing model details.
        
        
    Returns
    -------
    modObj : CSV_model object.
        
        
    Calls
    -----
    IO.CSV_Model.CSV_model
    
    """
    from toolkit.IO.netcdf_Model import netcdf_model
    
    wd = model.attrib["path"]
    model_fn = model.attrib["model_file_name"]
    x_var_name = model.attrib["x_var_name"]
    y_var_name = model.attrib["y_var_name"]
    values_var_name = model.attrib["values_var_name"]
    if "z_var_name" in model.attrib:
        z_var_name = model.attrib["z_var_name"]
    else:
        z_var_name = None
    #end if
    if "nanval" in model.attrib:
        nanval = fn.changetype(model.attrib["nanval"], float)
    else:
        nanval = np.nan
    #end if
    if "epsg" in model.attrib:
        proj_str = None
        epsg = fn.changetype(model.attrib["epsg"], int)
    else:
        proj_str = model.attrib["proj_str"]
        epsg = 0
    #end if
    if "swap_xy" in model.attrib:
        swap_xy = fn.changetype(model.attrib["swap_xy"], bool)
    else:
        swap_xy = False
    #end if
    if "interval" in model.attrib:
        interval = fn.changetype(model.attrib["interval"], int)
    else:
        interval = 1
    #end if
    
    modObj = netcdf_model(os.path.join(wd, model_fn), x_var_name, y_var_name, 
                          values_var_name, z_var_name=z_var_name, epsg=epsg, 
                          proj_str=proj_str, nanval=nanval, swap_xy=swap_xy,
                          interval=interval)
    
    modObj.epsg = epsg
    modObj.proj_str = proj_str
    return modObj
#end func

def read_geotiff_model(model):
    from toolkit.IO.geotiff_Model import geotiff_model

    wd = model.attrib["path"]
    model_fn = os.path.join(wd, model.attrib["model_file_name"])

    if "interval" in model.attrib:
        interval = fn.changetype(model.attrib["interval"], int)
    else:
        interval = 1
    #end if
    
    if 'epsg' in model.attrib:
        epsg = int(model.attrib['epsg'])
    else:
        epsg = False
    #end if
    modObj = geotiff_model(model_fn, interval=interval, force_epsg=epsg)

    return modObj
# end func
