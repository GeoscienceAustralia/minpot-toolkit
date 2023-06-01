"""
Deposit dataset input module
-----------------
This module is used to import deposit datasets from various input formats.
Supported formats are comma separated values (.csv), text (.txt), excel
spreadsheet (.xls, .xlsx).


Contacts
--------
Lachlan Adams
    - Lachlan.Adams@ga.gov.au or Lachlan.Adams.1996@outlook.com
    
    
"""

import os, time
import numpy as np
import pandas as pd
from geopandas.geoseries import GeoSeries
from shapely.geometry import Point
import toolkit.functions as fn
import geopandas

class DataManager():
    def __init__(self):
        """
        Class to facilitate the importation of mineral deposit datasets.
        
        
        """
        self.point_propertynames = []
        self.read_columns = []
        self.bounds = np.array([110, 155, -45, -10])
    #end func
    
    def ParseXMLNode(self, rootNode, commandList):
        """
        Function to get commands and the details of required datasets from xml 
        configuration file input.
        
        
        """
        for row in commandList:
            if hasattr(self, row[0]):
                tp = type(getattr(self, row[0]))
                val = fn.changetype(row[1], tp)
                setattr(self, row[0], val)
            #end if
        #end for
        self.datasets = rootNode.findall("deposits")
        self.newcolumns = rootNode.findall("newcolumn")
        self.properties = rootNode.findall("property")
        self.filter = rootNode.find("filter")
    #end func
    
    def Filter_Dataset(self, domain):
        """
        Function to filter input data to ensure it falls within the region of 
        coverage of the model being used.
        
        
        Calls
        -----
        filter_points
        
        
        """
        
        # Filter by location
        self.point_data_filt = filter_points(self.point_data, domain)
        if len(self.point_data_filt) == 0:
            raise Exception('No deposits found within region of coverage!')
    #end func
    
    def Make_Point_Property_Array(self):
        """
        Function to create an array which stores deposit attributes.
        
        
        Calls
        -----
        make_point_property_array
        
        
        """
        self.point_property_array = \
            make_point_property_array(self.point_data_filt, 
                                      self.point_propertynames)
    #end func
    
    def Read_Dataset(self, domain_dict):
        """
        Reads input files, filters input data, creates point property array, 
        and enters data into a dictionary to be used by other modules.
        
        
        """
        self.Read_Files()
        self.Filter_Dataset(domain_dict)
        self.Make_Point_Property_Array() 
        self.points_dict = {'point_data_filt': self.point_data_filt,
                            'point_property_array': self.point_property_array,
                            'n_deposits': len(self.point_data_filt)}
    #end func
    
    def Read_Files(self):
        """
        Reads input files and returns array of deposit locations and 
        properties.
        
        
        Calls
        -----
        read_deposit_dataset
        
        
        """
        self.point_data = read_deposit_dataset(self.datasets, 
                                               self.read_columns,
                                               self.point_propertynames,
                                               self.properties, self.filter, 
                                               self.newcolumns)
    #end func
#end class
    
class RandomPointsManager():
    """
    Class used to generate uniformly distributed random points, for analysis.
    
    
    """
    def __init__(self, **kwargs):
        """
        Initialises the object.
        
        
        """
        super().__init__(**kwargs)
        
        self.n_repeats_default = 100
        self.chunksize = 10000
        self.bounds = np.array([110, 155, -45, -10])
        
        # set attributes based on kwargs
        for key in kwargs.keys():
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
            #end if
        #end for
    #end func
    
    def ParseXMLNode(self, commandList):
        """
        Redefines default attributes based on input from xml configuration 
        file.
        
        
        """
        for row in commandList:
            if hasattr(self, row[0]):
                tp = type(getattr(self, row[0]))
                if tp == list:
                    val = fn.changetype(row[1], 'list', dtype='str')
                elif tp == np.ndarray:
                    val = fn.changetype(row[1], 'array', dtype='float')
                else:
                    val = fn.changetype(row[1], tp)
                #end if
                setattr(self, row[0], val)
            #end if
        #end for
    #end func
    
    def Prepare(self, domain_dict, points_dict):
        """
        Sets domain and n_deposits attributes for random point generation.
        
        
        """
        self.domain = domain_dict
        self.n_deposits = points_dict['n_deposits']
    #end func
    
    def seed_random_points(self, n_repeats=0):
        """
        Seeds uniformly distributed random points.
        

        Returns
        -------
        points_random : array
            Array of (lon, lat) locations of random points in equal earth 
            projection, with shape (n, 2)
            
            
        Calls
        -----
        filter_points
        
        functions.epsg_project
        
        
        """
        if n_repeats == 0:
            n_repeats = self.n_repeats_default
        #end if
        
        n_points = n_repeats*self.n_deposits
        
        # initialise new array to save points to
        random_points = np.empty((0, 2))
        
        while len(random_points) < n_points:
            random_points_chunk = \
                initialise_random_points_chunk(self.bounds, self.chunksize)
            points = filter_points(random_points_chunk, self.domain)
            random_points = np.vstack([random_points, points])
        #end while
        
        random_points = random_points[:n_points, :]
        
        random_points.resize((n_repeats, self.n_deposits, 2), refcheck=False)
        return random_points
    #end func
#end class
    
def filter_points(points, domain):
    """
    Filter points to be within a certain buffer of a station and within/outside
    polygons as defined by polygon_dict.
    

    Parameters
    ----------
    points : array
        Numpy array of x and y coordinates of points to be checked, with shape 
        (n, 2).
    
    domain : dictionary of shapely multipolygon
        Multipolygon describing the region of coverage for the model.
        

    Returns
    -------
    points : array
        Input points array filtered to contain only points within area of 
        interest and on land.
        
                       
    """
    
    coverage_grid = domain['grid']
    pointlist = GeoSeries([Point(pp) for pp in points[:, :2]])
    indices = np.arange(len(pointlist))
    
    print('Checking if', len(indices), 'points are within coverage grid')
    t0 = time.time()
    dist1 = pointlist.distance(coverage_grid)
    filt1 = dist1 <= 0.1
    indices = indices[filt1]
    print(len(indices), 'points, time =', time.time() - t0)
    
    return points[indices]
#end func
    
def initialise_random_points_chunk(bounds, chunksize):
    """
    Initialise points within the bounds as stored within the attribute.
    To ensure equal area representation across the globe they are seeded in 
    eastings and northings, then projected to longitude and latitude.
    

    Parameters
    ----------
    bounds : list of float
        minimum and maximum longitude and latitude values.
        
    chunksize : integer
        Number of points to be created per chunk.
        
        
    Returns
    -------
    random_locations_chunk : array
        Array of random location coordinates with shape (chunksize, 2).
        
    
    Calls
    -----
    functions.epsg_project
    
        
    """
    
    lon0, lon1, lat0, lat1 = bounds
    
    colat0 = np.pi*(90 - lat0)/180
    colat1 = np.pi*(90 - lat1)/180
    
    y0 = np.cos(colat0)
    y1 = np.cos(colat1)
    
    rand = np.random.random(size=chunksize)*(y1 - y0) + y0
    randcolat = np.arccos(rand)
    randlat = 90 - 180*randcolat/np.pi 
    randlon = np.random.random(size=chunksize)*(lon1 - lon0) + lon0
    
    return np.vstack([randlon, randlat]).T
#end func
    
def make_point_property_array(point_data_filt, point_propertynames):
    """
    Make and populate an array to contain point properties (filtered as for x, 
    y locations).
    
    
    Parameters
    ----------
    point_data_filt : array
        Array of coordinates of deposit locations for which distances to 
        contours are to be computed.
        
    point_propertynames : list of str
        Names of properties of the points used (e.g. deposit type).
    
    
    Returns
    -------
    point_property_array : array
        Structured array of point properties.
    
        
    """
    
    arrdtype = [(prop, np.float) for prop in point_propertynames]
    point_property_array = np.array(np.zeros(len(point_data_filt)),
                                         dtype=arrdtype)        
    for i, pname in enumerate(point_propertynames):
        point_property_array[pname] = point_data_filt[:, i+2]
    #end for
    return point_property_array
#end func
    
def read_csv(file, dataset, point_properties):
    """
    Function to read deposit information from a text file.
    
    
    Parameters
    ----------
    file : string
        Directory/name of input file.
        
    dataset : xml node
        ElementTree xml node containing information about the input dataset.
        
    point_properties : list of string
        Names of required columns for output dataframe. If column is not 
        present in input file, a column of zeros is created.
        
        
    Returns
    -------
    df : pandas dataframe
        DataFrame containing the coordinates and required point properties for
        each deposit in the input file.
        
    
    """
    columns = fn.changetype(dataset.attrib["columns"], list, dtype=str)
    df = pd.read_csv(file, names=columns).fillna(0)
    for prop in point_properties:
        if prop not in df.columns:
            df[prop] = np.zeros(len(df))
        #end if
    #end for
    return df[['longitude', 'latitude'] + point_properties]
#end func

def read_deposit_dataset(datasets, read_columns, point_properties, 
                         properties, rowfilter, newcolumns):
    """
    Reads the locations of mineral deposits from text files.
    
    
    Parameters
    ----------
    datasets : xml node
        ElementTree node of input parameters for the data sets.
        
    read_columns : list of string
        List of column names to be read from excel or csv files.
        
    point_properties : list of string
        List of column names for the array returned by the function. If the 
        point property is not present in the input data file, a column of zeros
        is created instead.
        
    properties : xml node
        ElementTree node of properties used to create a truth value array. For 
        example, a property may be "Cu (%) > 0".
        
    rowfilter : string
        Filter dataset by union or intersection of properties.
        
    newcolumns : xml node
        ElementTree node containing information to generate new columns for the
        data. For example, if we have Cu (%) and Tonnage (Mt) as columns, we 
        may make a new column "Cu (Mt) = Cu (%) * Tonnage (Mt)".
        
        
    Returns
    -------
    point_data : array
        Numpy array with shape (n, m+2) of coordinates for deposit locations 
        with the first and second coordinates being longitude and latitude 
        respectively, and the next m coordinates representing point properties.
        
        
    Calls
    -----
    read_txt
    
    read_csv
    
    read_excel
    
    add_columns
    
    filter_by_properties
    
    IO.XML.GetAttributeString
    
    
    """
    pointslist = list()
    for dataset in datasets:
        path = dataset.attrib["path"]
        file_name = dataset.attrib["file_name"]
        file = os.path.join(path, file_name)
        filetype = file_name.split('.')[-1]
        if filetype == 'txt':
            df = read_txt(file, dataset, read_columns)
        elif filetype == 'csv':
            df = read_csv(file, dataset, read_columns)
        elif filetype in ['xls', 'xlsx']:
            df = read_excel(file, dataset, read_columns)
        elif filetype in ['shp']:
            df = read_shp(file, dataset, read_columns)
        #end if
        pointslist.append(df)
    #end for
    point_df = pd.concat(pointslist, axis=0)
    point_df_newcols = add_columns(point_df, newcolumns)
    point_df_filt = filter_by_properties(point_df_newcols, properties, 
                                         rowfilter)
    return np.array(point_df_filt[['longitude', 'latitude'] + \
                                  point_properties])
#end func

def read_excel(file, dataset, point_properties):
    """
    Function to read deposit information from a text file.
    
    
    Parameters
    ----------
    file : string
        Directory/name of input file.
        
    dataset : xml node
        ElementTree xml node containing information about the input dataset.
        
    point_properties : list of string
        List of properties required for analysis.
        
        
    Returns
    -------
    df : pandas dataframe
        DataFrame containing the coordinates and required point properties for
        each deposit in the input file.
        
    
    """
    sheet_name = dataset.attrib["sheet_name"]
    lon_col_name = dataset.attrib["lon_col_name"]
    lat_col_name = dataset.attrib["lat_col_name"]
    df = pd.read_excel(file, sheet_name=sheet_name).fillna(0)
    for prop in point_properties:
        if prop not in df.columns:
            df[prop] = np.zeros(len(df))
        #end if
    #end for
    df = df[[lon_col_name, lat_col_name] + point_properties]
    df = df.rename(columns = {lon_col_name: 'longitude', 
                              lat_col_name: 'latitude'})
    return df
#end func
    
def read_txt(file, dataset, point_properties):
    """
    Function to read deposit information from a text file.
    
    
    Parameters
    ----------
    file : string
        Directory/name of input file.
        
    dataset : xml node
        ElementTree xml node containing information about the input dataset.
        
    point_properties : list of string
        Names of required columns for output dataframe. If column is not 
        present in input file, a column of zeros is created.
        
        
    Returns
    -------
    df : pandas dataframe
        DataFrame containing the coordinates and required point properties for
        each deposit in the input file.
        
    
    """
    columns = fn.changetype(dataset.attrib["columns"], list, dtype=str)
    points = np.genfromtxt(file, usecols=tuple(np.arange(len(columns))))
    df = pd.DataFrame(points, columns=columns).fillna(0)
    for prop in point_properties:
        if prop not in df.columns:
            df[prop] = np.zeros(len(df))
        #end if
    #end for
    return df[['longitude', 'latitude'] + point_properties]
#end func

def read_shp(file, dataset, point_properties):
    lon_col_name = dataset.attrib["lon_col_name"]
    lat_col_name = dataset.attrib["lat_col_name"]
    df = geopandas.read_file(file).fillna(0)
    for prop in point_properties:
        if prop not in df.columns:
            df[prop] = np.zeros(len(df))
        #end if
    #end for
    df = df[[lon_col_name, lat_col_name] + point_properties]
    df = df.rename(columns = {lon_col_name: 'longitude',
                              lat_col_name: 'latitude'})
    return df
# end func

def filter_by_properties(point_data, properties, rowfilter):
    """
    Calls
    -----
    functions.check_filter
    
    functions.check_property
    
    functions.get_indices
    
    
    """
    dct = {}
    keys = []
    for prop in properties:
        name = prop.attrib["name"]
        description = fn.changetype(prop.attrib["description"], list, 
                                    dtype=str)
        if description[0] in keys:
            prop1, condition, prop2 = description
            arr1 = dct[prop1]
            arr2 = dct[prop2]
            dct[name] = fn.check_filter(arr1, condition, arr2)
            keys.append(name)
        else:
            column, condition, value = description
            arr = np.array(point_data[column])
            try: 
                value = float(value)
            except:
                pass
            #end try
            dct[name] = fn.check_property(arr, condition, value)
            keys.append(name)
        #end if
    #end for
    keys = fn.changetype(rowfilter.attrib["filters"], list, dtype=str)
    condition = rowfilter.attrib["condition"]
    arrlist = [dct[key] for key in keys]
    indices = fn.get_indices(arrlist, condition)
    return point_data[indices]
#end func
    
def add_columns(point_data, newcolumns):
    """
    Calls
    -----
    column_operation
    
    
    """
    for newcolumn in newcolumns:
        name = newcolumn.attrib["name"]
        column1, operation, value = \
            fn.changetype(newcolumn.attrib["description"], list, dtype=str)
        column1 = point_data[column1]
        if value in point_data.columns:
            column2 = point_data[value]
            point_data[name] = column_operation(column1, column2, operation)
        else:
            point_data[name] = column_operation(column1, float(value), 
                                                operation)
        #end if
    return point_data
#end func
    
def column_operation(x, y, operation):
    """
    Function to operate on two columns of a dataframe based on xml input 
    command.
    
    
    """
    if operation == '*':
        return x*y
    elif operation == '+':
        return x + y
    elif operation == '-':
        return x - y
    elif operation == '/':
        return x/y
    #end if
#end func