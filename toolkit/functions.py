"""
Functions
-----------------
This module contains a set of functions used throughout the mineral potential
toolkit.


Contacts
--------
Lachlan Adams
    - Lachlan.Adams@ga.gov.au or Lachlan.Adams.1996@outlook.com
    
    
"""

import numpy as np
from scipy.special import factorial
from pyproj import Transformer

def changetype(var, tp, dtype="str", sep=', '):
    """
    Function to convert a variable to a different type. Tested for input
    strings only.
    
    
    Parameters
    ----------
    var : string, float, integer, or boolean
        Variable to convert to a different type.
        
    tp : string
        Name of type to convert variable to.
        
    dtype : string
        Data type for output if type is array or list.
        
    sep : string
        Separator for if the variable is to be converted into an array or list.
        
        
    Returns
    -------
    var : string, float, integer, boolean, array, or list
        Converted variable.
        
        
    """
    if tp == 'str' or tp == str:
        return str(var)
    elif tp == 'float' or tp == float:
        return float(var)
    elif tp == 'int' or tp == int:
        return int(var)
    elif tp == 'bool' or tp == bool:
        return var == 'True' or var == '1'
    elif tp == 'list' or tp == list:
        if dtype == "str" or dtype == str:
            return var.split(sep)
        elif dtype == "float" or dtype == float:
            return [float(item) for item in var.split(sep)]
        elif dtype == "int" or dtype == int:
            return [int(item) for item in var.split(sep)]
        elif dtype == "bool" or dtype == bool:
            return [(item == "True" or item == '1') for item in var.split(sep)]
        #end if
    elif tp == 'array' or tp == np.ndarray:
        if dtype == "str" or dtype == str:
            dtype = "float"
        #end if
        if dtype == "float" or dtype == float:
            return np.array([float(item) for item in var.split(sep)])
        elif dtype == "int" or dtype == int:
            return np.array([int(item) for item in var.split(sep)])
        elif dtype == "bool" or dtype == bool:
            return np.array([(item == "True" or item == '1') \
                             for item in var.split(sep)])
        #end if
    else:
        return None
    #end if
#end func

def check_filter(arr1, condition, arr2):
    """
    Function to determine if 'arr1' and 'arr2' fulfil 'condition' at every 
    coordinate.
    
    
    Parameters
    ----------
    arr1 : numpy array
        Boolean array
        
    arr2 : numpy array
        Boolean array
        
    condition : string
        Logical operator used to check arr1 and arr2. Currently supported are
        'or', 'and', 'xor', and 'nand'.
        
    
    Returns
    -------
    arr : numpy array
        Boolean array.
    
    
    """
    if condition == 'or':
        return np.logical_or(arr1, arr2)
    elif condition == 'and':
        return np.logical_and(arr1, arr2)
    elif condition == 'xor':
        return np.logical_xor(arr1, arr2)
    elif condition == 'nand':
        return ~np.logical_and(arr1, arr2)
    #end if
#end func
    
def check_property(arr, condition, value):
    """
    Function used to compare 'arr' to 'value'.
    
    
    Parameters
    ----------
    arr : numpy array
        Array of values to compare to 'value'.
        
    value : float
        Value to compare 'arr' to.
        
    condition : string
        Condition used to compare 'arr' and 'value'. Currently supported are
        '>', '<', '==', '!=', '>=', '<='.
        
    
    Returns
    -------
    arr1 : numpy array
        Boolean array representing whether the corresponding coordinate in 
        'arr' fulfilled 'condition' when compared to 'value'.
        
    
    """
    if condition == '>':
        return arr > value
    elif condition == '<':
        return arr < value
    elif condition == '==':
        return arr == value
    elif condition == '!=':
        return arr != value
    elif condition == '>=':
        return arr >= value
    elif condition == '<=':
        return arr <= value
    #end if
#end func
    
def epsg_project(x, y, epsg_from, epsg_to, proj_str=None):
    """
    Project (x, y) points using the pyproj modules.
    

    Parameters
    ----------
    x : float
        x coordinate of point.
   
    y : float
        y coordinate of point.
    
    epsg_from : int
        EPSG code of (x, y) points provided. To provide custom projection, set 
        to 0 and provide proj_str.
    
    epsg_to : int
        EPSG code to project to. To provide custom projection, set to 0 and 
        provide proj_str.
    
    
    Optional parameters
    -------------------
    proj_str : str
        Proj4 string to provide to pyproj if using custom projection. This proj 
        string will be applied if epsg_from or epsg_to == 0. The default is 
        None.
        

    Returns
    -------
    xp : float
        x coordinate of projected point.
        
    yp : float
        y coordinate of projected point.
        
    
    """
    
    proj_list = []
    for epsg in [epsg_from, epsg_to]:
        if epsg == 0:
            #proj_list.append(Proj(proj_str))
            proj_list.append(proj_str)
        else:
            #proj_list.append(pyproj.Proj(init='epsg:%1i'%epsg))
            #proj_list.append(Proj('epsg:%1i'%epsg))
            proj_list.append('epsg:%1i'%epsg)
        #end if
    #end for

    p1, p2 = proj_list
    t = Transformer.from_crs(p1, p2, always_xy=True)

    #return transform(p1, p2, x, y, always_xy=True)
    return t.transform(x, y)
#end func
    
def get_indices(properties, condition):
    if condition == 'union':
        return np.any(properties, axis=0)
    elif condition == 'intersection':
        return np.all(properties, axis=0)
    #end if
#end func
    
def gradient(x, y, z, values):
    """
    Function to return the gradient of a slice of a rank-3 tensor with a grid
    of values at coordinates x, y, z.
    
    
    Parameters
    ----------
    values : array
        Rank-3 tensor of function values.
    x : array
        Matrix of x coordinates for values grid.
    y : array
        Matrix of y coordinates for values grid.
    z : array
        Vector of z coordinates for values grid.
        
        
    Returns
    -------
    gradient_values : array
        Rank-3 tensor of gradient values.
        
        
    """
    gradient_values = np.zeros_like(values)
    gx = np.zeros_like(x)
    gy = np.zeros_like(y)
    
    for i in range(len(z)):
        # Gradient in y direction
        # Bottom edge (forward difference method)
        gy[0,:] = (values[1,:,i] - values[0,:,i])/(y[1,:] - y[0,:])
        
        # Top edge (backward difference method)
        gy[-1,:] = (values[-1,:,i] - values[-2,:,i])/(y[-1,:] - y[-2,:])
        
        # Middle (central difference method)
        gy[1:-1,:] = (values[2:,:,i] - values[:-2,:,i])/(y[2:,:] - y[:-2,:])
        
        # Gradient in x direction
        # Left edge (forward difference method)
        gx[:,0] = (values[:,1,i] - values[:,0,i])/(x[:,1] - x[:,0])
        
        # Right edge (backward difference method)
        gx[:,-1] = (values[:,-1,i] - values[:,-2,i])/(x[:,-1] - x[:,-2])
        
        # Middle (central difference method)
        gx[:,1:-1] = (values[:,2:,i] - values[:,:-2,i])/(x[:,2:] - x[:,:-2])
        
        # Scalar
        g = np.sqrt(gx**2 + gy**2)
        ind = np.logical_and(np.isfinite(values[:,:,i]), ~np.isfinite(g))
        g[ind] = 0.0
        gradient_values[:,:,i] = g
    #end for
    return gradient_values
#end func
    
def gradient_sph(lon, lat, r, values):
    """
    Function to return the gradient of a slice of a rank-3 tensor with a grid
    of values in spherical coordinates lon, lat, r.
    
    
    Parameters
    ----------
    values : array
        Rank-3 tensor of function values.
    lon : array
        Matrix of longitude coordinates for values grid.
    lat : array
        Matrix of latitude coordinates for values grid.
    r : array
        Vector of radii for values grid.
        
        
    Returns
    -------
    gradient_values : array
        Rank-3 tensor of gradient values.
        
        
    """
    gradient_values = np.zeros_like(values)
    degrad = np.pi/180
    theta = (90 - lat)*degrad
    phi = lon*degrad
    df_dtheta = np.zeros_like(theta)
    df_dphi = np.zeros_like(phi)
    
    for i in range(len(r)):
        # Gradient in theta direction
        # Bottom edge (forward difference method)
        df_dtheta[0,:] = (values[1,:,i] - values[0,:,i])/ \
                         (theta[1,:] - theta[0,:])
        
        # Top edge (backward difference method)
        df_dtheta[-1,:] = (values[-1,:,i] - values[-2,:,i])/ \
                       (theta[-1,:] - theta[-2,:])
        
        # Middle (central difference method)
        df_dtheta[1:-1,:] = (values[2:,:,i] - values[:-2,:,i])/ \
                       (theta[2:,:] - theta[:-2,:])
        
        # Gradient in phi direction
        # Left edge (forward difference method)
        df_dphi[:,0] = (values[:,1,i] - values[:,0,i])/(phi[:,1] - phi[:,0])
        
        # Right edge (backward difference method)
        df_dphi[:,-1] = (values[:,-1,i] - values[:,-2,i])/ \
                        (phi[:,-1] - phi[:,-2])
        
        # Middle (central difference method)
        df_dphi[:,1:-1] = (values[:,2:,i] - values[:,:-2,i])/ \
                       (phi[:,2:] - phi[:,:-2])
        
        # Scalar
        gtheta = df_dtheta/r[i]
        s = np.sin(theta)
        s[s == np.sin(np.pi)] = 0.0
        gphi = df_dphi/(r[i]*s)
        g = np.sqrt(gtheta**2 + gphi**2)
        ind = np.logical_and(np.isfinite(values[:,:,i]), ~np.isfinite(g))
        g[ind] = 0.0
        gradient_values[:,:,i] = g
    #end for
    return gradient_values
#end func
    
def nearest_index(val, array):
    """
    Find the index of the nearest value in an array to input value.
    
    
    Parameters
    ----------
    val : float
        The value to search for.
    
    array : array
        The array to search in.
        
    
    Returns
    -------
    index : int
        Index of nearest value in array.
        
    
    """
    # absolute difference between value and array
    diff = np.abs(array - val)
    
    return np.where(diff==min(diff))[0][0]
#end func
    
def normal(x, m, s):
    """
    Computes the normal probability density value at x given mean m and 
    standard deviation s.
    
    
    """
    return np.sqrt(np.exp(-((x - m)/s)**2)/(2*np.pi))/s
#end func
    
def Poisson(x, m):
    """
    Returns the value of the Poisson probability density function at x given
    mean m.
    
    
    """
    return m**x * np.exp(-m)/factorial(x)
#end func
    
def binomial(x, m, n):
    """
    Returns the value of the binomial probability density function at x given
    mean m and number of points n.
    
    
    """
    return factorial(n)*(m/n)**x*(1-(m/n))**(n-x)/(factorial(n-x)*factorial(x))
#end for