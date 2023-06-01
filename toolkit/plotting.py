"""
Plotting module
-----------------
This module is used to generate plots representing the results of tests 
performed using the statistics module.


Contacts
--------
Lachlan Adams
    - Lachlan.Adams@ga.gov.au or Lachlan.Adams.1996@outlook.com
    
    
"""

import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import toolkit.functions as fn

from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from shapely.geometry import Point, Polygon, LineString

def plot_distribution_at_depth_slice(obj, depths, wd_images):
    """
    Generates plots showing the expected distribution of random locations, and 
    the distribution of real deposits, at each depth slice when performing the
    Poisson or Binomial p test.
    
    
    Parameters
    ----------
    obj: stat_corr.(various) object
        Object containing the results of a statistical test.
        
    depths: numpy array
        1D array of depth values for distribution to be plotted at.
        
    wd_images: string
        Directory in which to save figure.
    
    
    Calls
    -----
    plot_distribution
    
    
    """
    path = os.path.join(wd_images, obj.name)
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    if type(depths) == float or type(depths) == int:
        plot_distribution(obj, depths)
        filename = \
            os.path.join(path, '%s_km_depth_distribution.png'%(depths/1000))
        plt.savefig(filename)
        plt.show()
    elif type(depths) == list or type(depths) == np.ndarray:
        for depth in depths:
            plot_distribution(obj, depth)
            filename = \
                os.path.join(path, '%s_km_depth_distribution.png'%(depth/1000))
            plt.savefig(filename)
            plt.show()
        #end for
    #end if
    return
#end func
    
def plot_distribution(obj, depth):
    """
    Generates a plot showing a histogram of the number of random locations 
    which fell within an area of interest compared to those expected if the
    Poisson or Binomial distribution is followed.
    
    
    Parameters
    ----------
    obj: stat_corr.(various) object
        Object containing the results of a statistical test.
        
    depth: float
        Depth at which to plot the distribution.
    
    
    Calls
    -----    
    as_si
    
    functions.nearest_index
    
    functions.Poisson
    
    functions.Binomial
    
    
    """
    didx = fn.nearest_index(depth, obj.distance_arrays[0]['depth'])
    n_points = obj.n_points
    rvals = obj.rvals[didx]*n_points
    n_repeats = len(rvals)
    x = obj.xvals[didx]*n_points
    m = np.mean(rvals)
    maxbin = int(np.min([n_points+1, np.max([x*2, np.max(rvals)*2])]))
    freqs, binedges = np.histogram(rvals, bins=maxbin, range=(0, maxbin))
    xvals = binedges[:-1]
    plt.figure()
    plt.plot(xvals, freqs, c='b', label='Histogram')
    if obj.distribution == 'Poisson':
        textstr = \
            'P(x = %s) = %s'%(x, "\n${0:s}$".format(as_si(fn.Poisson(x, m))))
        plt.plot(xvals, fn.Poisson(xvals, m)*n_repeats, c='r', label='Poisson')
    elif obj.distribution == 'Binomial':
        textstr = 'P(x = %s) = %s'%(x, "\n${0:s}$".format(as_si(\
                    fn.binomial(x, m, n_points))))
        plt.plot(xvals, fn.binomial(xvals, m, n_points)*n_repeats, c='r', 
                 label='Binomial')
    #end if
    plt.axvspan(x-0.1, x+0.1, ymin=0, ymax=0.25, color='green', alpha=0.5, 
                label='Observed')
    plt.text(x+0.2, 0.2*np.max(freqs), textstr, ha='left', 
             va='center', color='k', fontsize=10) 
    plt.legend()
    plt.xlabel('Number of points in region')
    tick_int = int(np.ceil(maxbin/10))
    plt.xticks(binedges[::tick_int])
    plt.ylabel('Frequency')
    plt.title('%s km Distribution function'%(depth/1000))
    return
#end func
    
def plot_Poisson_p_values(obj, path):
    """
    Generates plots showing the differences between distributions for random 
    locations versus deposit points, and plots showing the p values for the 
    Poisson p test at each depth.
    
    
    Parameters
    ----------
    obj: stat_corr.(various) object
        Object containing the results of a statistical test.
        
    path: string
        Directory in which to save figure.
    
    
    """
    depths = np.array(obj.distance_arrays[0]['depth'])
    rvals = obj.rvals
    xvals = obj.xvals
    mvals = np.mean(rvals, axis=1)
    svals = np.std(rvals, axis=1)
    minvals = np.max([np.zeros_like(mvals), mvals - svals], axis=0)
    maxvals = mvals + svals
    filename = os.path.join(path, 'Points_within_region.png')
    plt.figure()
    plt.plot(depths/1e3, xvals*100, label='real')
    plt.fill_between(depths/1e3, minvals*100, maxvals*100, alpha=0.5, 
                     label='random')
    plt.xlabel('Depth (km)')
    plt.ylabel('Points within region (%)')
    plt.title('Real vs random results')
    plt.legend()
    plt.savefig(filename)
    plt.show()
    
    filename = os.path.join(path, 'P_values.png')
    plt.figure()
    plt.plot(depths/1e3, obj.pvals, c='b', label='Poisson p values')
    plt.plot(depths/1e3, 1-obj.chi_p_vals, c='r', label='\chi^2 p values')
    plt.yscale('log')
    plt.xlabel('Depth (km)')
    plt.ylabel('p value')
    plt.ylim([None, 1])
    plt.title('Poisson p test results')
    plt.legend()
    plt.savefig(filename)
    plt.show()
    return
#end func
    
def plot_Binomial_p_values(obj, path):
    """
    Generates plots showing the differences between distributions for random 
    locations versus deposit points, and plots showing the p values for the 
    Poisson p test at each depth.
    
    
    Parameters
    ----------
    obj: stat_corr.(various) object
        Object containing the results of a statistical test.
        
    path: string
        Directory in which to save figure.
        
    
    """
    depths = np.array(obj.distance_arrays[0]['depth'])
    rvals = obj.rvals
    xvals = obj.xvals
    mvals = np.mean(rvals, axis=1)
    svals = np.std(rvals, axis=1)
    minvals = np.max([np.zeros_like(mvals), mvals - svals], axis=0)
    maxvals = mvals + svals
    filename = os.path.join(path, 'Points_within_region.png')
    plt.figure()
    plt.plot(depths/1e3, xvals*100, label='real')
    plt.fill_between(depths/1e3, minvals*100, maxvals*100, alpha=0.5, 
                     label='random')
    plt.xlabel('Depth (km)')
    plt.ylabel('Points within region (%)')
    plt.title('Real vs random results')
    plt.legend()
    plt.savefig(filename)
    plt.show()
    
    filename = os.path.join(path, 'P_values.png')
    plt.figure()
    plt.plot(depths/1e3, obj.pvals, c='b', label='Binomial p values')
    plt.plot(depths/1e3, 1-obj.chi_p_vals, c='r', label='\chi^2 p values')
    plt.yscale('log')
    plt.xlabel('Depth (km)')
    plt.ylabel('p value')
    plt.ylim([None, 1])
    plt.title('Binomial p test results')
    plt.legend()
    plt.savefig(filename)
    plt.show()
    return
#end func

def as_si(x, ndp=0):
    """
    Formats a number in scientific notation i.e. x 10 -3 etc.
    Modified from https://stackoverflow.com/questions/31453422/displaying-numbers-with-x-instead-of-e-scientific-notation-in-matplotlib/31453961.
    

    Parameters
    ----------
    x : float
        Number to be represented.
        
    ndp: integer
        Number of decimal points in representation.
        

    Returns
    -------
    x : string
        Input number represented in scientific notation
        
    
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times10^{{{e:d}}}'.format(m=m, e=int(e))
#end func
    
def plot_cdf(obj, depth, normed=True, colour='b', pad=20e3):
    
    """
    Function to plot the cumulative distribution function of distances to a 
    target contour from deposit locations at a given target depth.
    
    
    Parameters
    ----------
    obj : statistical_correlation.cdf object
    
    depth : float
        Target depth for the CDF to be plotted at.
        
    normed : boolean
        Flag to plot (or not plot) normed CDF as opposed to CDF.
        
    colour : string
        Colour of plot.
        
    pad : float
        Left and right limits (metres) of plot outside of distance range .
        
    
    Calls
    -----
    functions.nearest_index
    
    
    """
    
    didx = fn.nearest_index(depth, obj.cdf_array['depth'])
        
    # get data to plot
    if normed:
        cdf, cdfr, cdfstd = obj.cdf_normed[didx]*100, \
            obj.cdf_random_normed[didx]*100, \
            obj.cdf_random_std[didx]*100
        idx1 = np.where(cdfr==100)[0][0]
        y0, y1 = -5, 105
    else:
        cdf, cdfr, cdfstd = obj.cdf_array['cdf'][didx],\
            np.mean(obj.cdf_array['cdf_random'][didx], axis=0),\
            np.std(obj.cdf_array['cdf_random'][didx], axis=0)
        y0, y1 = np.array([-0.05, 1.05])*obj.n_deposits
        idx1 = np.where(cdfr==obj.n_deposits)[0][0]
    #end if
    
    idx0 = np.where(cdfr==0)[0][-1]
    
    x0, x1 = obj.cdf_array['binc'][didx,idx0] - pad, \
             obj.cdf_array['binc'][didx,idx1] + pad
        
    plt.plot(obj.cdf_array['binc'][didx]/1000, cdf, color=colour, lw=2.0, label='Raw deposits')
    
    plt.fill_between(obj.cdf_array['binc'][didx]/1000,
                     np.max([np.zeros_like(cdfr), cdfr - cdfstd], axis=0), 
                     cdfr + cdfstd,
                     color=colour,
                     lw=0.0,
                     alpha=0.5,
                     label='Random locations')
    plt.plot(obj.cdf_array['binc'][didx]/1000, 
             cdfr, 
             color=colour,
             lw=0.5)
    plt.xlim(x0/1000, x1/1000)
    
    plt.legend()
    plt.title('%s km depth'%(depth/1000))
    plt.xlabel('Distance to contour (km)')
    plt.ylabel('CDF (%)')
    return
#end func

def plot_cdf_at_depth_slice(depths, wd_images, normed=True, dct={}, 
                            feature=None, obj=None):
    
    """
    Function to plot the cumulative distribution function of distances to a 
    target contour from deposit locations at a range of target depths.
    
    
    Parameters
    ----------
    obj : stat_corr.Kolmogorov_Smirnov_p_test object
    
    dct : dictionary
        Dictionary containing stat_corr.Kolmogorov_Smirnov_p_test 
        objects for which plots are needed. Use if obj is not defined.
        
    feature : string
        Name of feature to plot CDF for. Use if obj is not defined.
        
    depths : array
        Target depths for plotting.
        
    wd_images : string
        Directory for image to be saved in.
        
        
    Calls
    -----
    plot_cdf
    
    
    """
    if obj is None:
        obj = dct[feature]
        path = os.path.join(wd_images, feature)
    else:
        path = os.path.join(wd_images, obj.name)
    #end if
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    if type(depths) == float or type(depths) == int:
        plot_cdf(obj, depths, normed=normed)
        filename = os.path.join(path, '%s_km_depth_CDF.png'%(depths/1000))
        plt.savefig(filename)
        plt.show()
    elif type(depths) == list or type(depths) == np.ndarray:
        for depth in depths:
            plot_cdf(obj, depth, normed=normed)
            filename = os.path.join(path, '%s_km_depth_CDF.png'%(depth/1000))
            plt.savefig(filename)
            plt.show()
        #end for
    #end if
    return
#end func

def plot_heat_map(wd_images, ks_depths=['min'], dct={}, feature=None, 
                  obj=None):
    """
    Function to create a "heat map" of "D" values at each depth slice for a
    CDF object. "D" is the distance between real cumulative distribution
    function and mean cumulative distribution for random locations.
    
    
    Parameters
    ----------
    obj : stat_corr.Kolmogorov_Smirnov_p_test object.
    
    dct : dictionary
        Dictionary containing stat_corr.Kolmogorov_Smirnov_p_test 
        objects for which plots are needed. Use if obj is not defined.
        
    feature : string
        Name of feature to plot heat map for. Use if obj is not defined.
    
    wd_images : string
        Directory for image to be saved in.
    
    ks_depths : list
        List of depths to show the Kolmogorov-Smirnov P value as text overlaid
        on top of the heat map.
        
                    
    Returns
    -------
    None
    
        
    Calls
    -----
    as_si
    
    functions.nearest_index
    
        
    """
    if obj is None:
        obj = dct[feature]
        path = os.path.join(wd_images, feature)
    else:
        path = os.path.join(wd_images, obj.name)
    #end if
    bins = obj.cdf_array['bins']/1000.0
    depths = obj.cdf_array['depth']/1000.0
    cdfr = np.mean(obj.cdf_array['cdf_random'][1], axis=0)
    idx0 = np.where(cdfr==0)[0][-1]
    idx1 = np.where(cdfr==np.max(cdfr))[0][0]
    dvalues = obj.dvalues
    ks = obj.kolmogorov_smirnov
    dmax = obj.dmax
    
    ax = plt.subplot(111)
    im=ax.pcolor(bins, depths, dvalues, vmin=0, vmax=1, cmap='hot_r')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Depth (km)')
    plt.colorbar(im, ax=ax)
    ax.set_xlim(bins[1][idx0], bins[1][idx1])
    ax.set_ylim(max(depths), 0)
    
    ax2 = plt.twiny()
    for ks_depth in ks_depths:
        if ks_depth == 'min':
            didx = np.where(ks == np.amin(ks))[0][0]
        else:
            didx = fn.nearest_index(ks_depth,depths)
        #end if
        pkslabel = r"P$_{KS}=$"+"\n${0:s}$".format(as_si(ks[didx]))
        ax2.text(dmax[didx], depths[didx], pkslabel,
                 ha='left', va='center', color='k', fontsize=10)   
    #end for

    ax2.plot(dmax, depths, 'k-')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(max(depths), 0)
    ax2.set_xlabel('D_max')
        
    plt.title('Heat plot')   
    
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    filename = os.path.join(path, 'Heatmap.png')
    plt.savefig(filename)
    plt.show()
    return
#end func

def plot_kolmogorov_smirnov_values(wd_images, dct={}, feature=None, obj=None, save=True):
    
    """
    Function to create a line plot of the Kolmogorov-Smirnov P value as a 
    function of depth.
    
    
    Parameters
    ----------
    obj : stat_corr.Kolmogorov_Smirnov_p_test object.
    
    dct : dictionary
        Dictionary containing stat_corr.Kolmogorov_Smirnov_p_test 
        objects for which plots are needed. Use if obj is not defined.
        
    feature : string
        Name of feature to make plot for. Use if obj is not defined.
    
    wd_images : string
        Directory for image to be saved in.
        
        
    """
    if obj is None:
        obj = dct[feature]
        path = os.path.join(wd_images, feature)
    else:
        path = os.path.join(wd_images, obj.name)
    #end if
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    filename = os.path.join(path, 'Kolmogorov_Smirnov_values.png')
    plt.plot(obj.cdf_array['depth']/1000, obj.kolmogorov_smirnov)
    plt.yscale('log')
    plt.xlabel('Depth (km)')
    plt.ylabel('p_{ks} value')
    plt.ylim([None, 1])
    plt.title('Kolmogorov-Smirnov test results')
    if save:
        plt.savefig(filename)
    plt.show()
    return
#end func

def plot_feature(feature, points_dict, depth, wd_images, bgimage="", save=True):
    """
    Function to create a representation of the modelled resistivity at a 
    specified depth below the surface, overlaid on a map. Also includes 
    locations of known mineral deposits as a scatter plot.
    
    
    Parameters
    ----------
    feature : feat_eng.feature_object object.
    
    points_dict : dictionary
        Dictionary which must contain 'point_data_filt', an array of deposit
        locations.
    
    depth : float
        Target depth value.
        
    wd_images : string
        Directory for image to be saved in.
        
    bgimage : string or matplotlib.imread object
        Background image on which to overlay the feature.
        
        
    Calls
    -----
    functions.nearest_index
    
    
    """
    gclon = feature.gclon
    gclat = feature.gclat
    gcz = feature.gcz
    depth_slice = fn.nearest_index(depth, gcz)
    feature_values = feature.feature_values
    point_data_filt = points_dict['point_data_filt']
    path = os.path.join(wd_images, feature.name)
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    filename = os.path.join(path, \
        '%s_km_%s_map.png'%(depth/1000, feature.name))
    lonmin = np.floor(np.min(gclon))
    lonmax = np.ceil(np.max(gclon))
    latmin = np.floor(np.min(gclat))
    latmax = np.ceil(np.max(gclat))
    plt.figure(figsize=(10,7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())
    if type(bgimage) != str:
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90], 
                   zorder=1)
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90], 
                   alpha=0.5, 
                   zorder=3)
    #end if
    vals = feature_values[:, :, depth_slice]
    vmin = -1
    vmax = 1
    plt.pcolormesh(gclon, gclat, vals, cmap='RdBu', vmin=vmin, vmax=vmax, 
                   zorder=2)
    points = plt.scatter(point_data_filt[:,0], point_data_filt[:,1], c='r', 
                         zorder=4, label='Deposit')
    #red_patch = mpatches.Patch(color='red', label='Not feature')
    blue_patch = mpatches.Patch(color='blue', label='Feature')
    plt.legend(handles=[points, blue_patch]) #, red_patch])
    ax.coastlines(zorder=5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    #plt.colorbar(im2, ax=ax, label='Feature T/F')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('%s km %s'%(depth/1000, feature.name))
    if save:
        plt.savefig(filename)
    #end if
    plt.show()
    return
#end func

def plot_feature_at_depth_slice(feature, points_dict, depths, wd_images, 
                                bgimage=None):
    """
    Function to plot a feature at a series of depth slices.
    
    
    Parameters
    ----------
    feature : feat_eng.feature_object object.
    
    points_dict : dictionary
        Dictionary which must contain 'point_data_filt', an array of deposit
        locations.
    
    depths : numpy array
        Target depth values.
        
    wd_images : string
        Directory for image to be saved in.
        
    bgimage : string or matplotlib.imread object
        Background image on which to overlay the feature.
    
    
    Calls
    -----
    plot_feature
    
    
    """
    if type(depths) in [float, int]:
        plot_feature(feature, points_dict, depths, wd_images, bgimage=bgimage)
    else:
        for depth in depths:
            plot_feature(feature, points_dict, depth, wd_images, 
                         bgimage=bgimage)
        #end for
    #end if
#end func
    
def plot_pdf(obj, depth, normed=True, colour='b', pad=20e3):
    
    """
    Function to plot the cumulative distribution function of distances to a 
    target contour from deposit locations at a given target depth.
    
    
    Parameters
    ----------
    obj : stat_corr.Kolmogorov_Smirnov_p_test object
    
    depth : float
        Target depth for the CDF to be plotted at.
        
    normed : boolean
        Flag to plot (or not plot) normed PDF as opposed to PDF.
        
    colour : string
        Colour of plot.
        
    pad : float
        Left and right limits (metres) of plot outside of distance range .
        
    
    Calls
    -----
    functions.nearest_index
    
    
    """
    
    didx = fn.nearest_index(depth, obj.cdf_array['depth'])
        
    # get data to plot
    if normed:
        pdf, pdfr, pdfstd = obj.pdf_normed[didx]*100, \
            obj.pdf_random_normed[didx]*100, \
            obj.pdf_random_std[didx]*100
        y0, y1 = -5, 105
    else:
        pdf, pdfr, pdfstd = obj.cdf_array['pdf'][didx],\
            np.mean(obj.cdf_array['pdf_random'][didx], axis=0),\
            np.std(obj.cdf_array['pdf_random'][didx], axis=0)
        y0, y1 = np.array([-0.05, 1.05])*obj.n_deposits
    #end if
    
    x0, x1 = np.min(obj.cdf_array['binc']) - pad, \
             np.max(obj.cdf_array['binc']) + pad
        
    plt.plot(obj.cdf_array['binc'][didx]/1000, pdf, color=colour)
    plt.fill_between(obj.cdf_array['binc'][didx]/1000,
                     np.max([pdfr - pdfstd, np.zeros_like(pdfr)], axis=0), 
                     pdfr + pdfstd,
                     color=colour,
                     alpha=0.5)
    plt.xlim(x0/1000, x1/1000)
    
    plt.title('%s km depth'%(depth/1000))
    plt.xlabel('Distance to contour (km)')
    plt.ylabel('PDF (%)')
    return
#end func

def plot_pdf_at_depth_slice(depths, wd_images, normed=True, dct={}, 
                            feature=None, obj=None):
    
    """
    Function to plot the cumulative distribution function of distances to a 
    target contour from deposit locations at a range of target depths.
    
    
    Parameters
    ----------
    obj : stat_corr.Kolmogorov_Smirnov_p_test object
    
    dct : dictionary
        Dictionary of stat_corr.Kolmogorov_Smirnov_p_test objects.
        Use if obj is not defined.
        
    feature : string
        Name of feature to plot PDF for. Use if obj is not defined.
        
    depths : array
        Target depths for plotting.
        
    wd_images : string
        Directory for image to be saved in.
        
        
    Calls
    -----
    plot_pdf
    
    
    """
    if obj is None:
        obj = dct[feature]
        path = os.path.join(wd_images, feature)
    else:
        path = os.path.join(wd_images, obj.name)
    #end if
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    if type(depths) == float or type(depths) == int:
        plot_pdf(obj, depths, normed=normed)
        filename = os.path.join(path, '%s_km_depth_PDF.png'%(depths/1000))
        plt.savefig(filename)
        plt.show()
    elif type(depths) == list or type(depths) == np.ndarray:
        for depth in depths:
            plot_pdf(obj, depth, normed=normed)
            filename = os.path.join(path, '%s_km_depth_PDF.png'%(depth/1000))
            plt.savefig(filename)
            plt.show()
        #end for
    #end if
    return
#end func
    
def plot_quantile_map(obj, bgimage="", wd_images='.'):
    """
    Generates a plot of the quantile which a grid point falls in based on its
    distance to a feature, in terms of the distribution of distances from 
    deposits to the feature.
    
    
    Parameters
    ----------
    obj : pros_mod.prospectivity_map object
    
    bgimage : string of matplotlib.imread object
        Background image on which to overlay plot.
        
    wd_images : string
        Directory in which to save image.
        
    
    """
    lon = obj.lon
    lat = obj.lat
    values = obj.quantile_values
    lonmin = np.floor(np.min(lon))
    lonmax = np.ceil(np.max(lon))
    latmin = np.floor(np.min(lat))
    latmax = np.ceil(np.max(lat))
    plt.figure(figsize=(10,7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())
    if type(bgimage) != str:
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90], 
                   zorder=1)
    #end if
    vmin = 0
    vmax = 1
    im2 = plt.pcolormesh(lon, lat, values[:,:,0], cmap='hot_r', vmin=vmin, 
                         vmax=vmax, zorder=2)
    ax.coastlines(zorder=4)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    plt.colorbar(im2, ax=ax, label='Quantile')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Quantile map')
    name = obj.stat_obj.name
    path = os.path.join(wd_images, name)
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    filename = os.path.join(path, 'Quantile_map.png')
    plt.savefig(filename)
    plt.show()
    return
#end func
    
def plot_pdf_map(obj, bgimage="", wd_images='.'):
    """
    Generates a plot of the 2D spatial PDF for a distribution of deposits in 
    relation to their distance from a feature.
    
    
    Parameters
    ----------
    obj : pros_mod.prospectivity_map object
    
    bgimage : string of matplotlib.imread object
        Background image on which to overlay plot.
        
    wd_images : string
        Directory in which to save image.
        
    
    """
    lon = obj.lon
    lat = obj.lat
    values = obj.pdf_values/np.nanmax(obj.pdf_values)
    lonmin = np.floor(np.min(lon))
    lonmax = np.ceil(np.max(lon))
    latmin = np.floor(np.min(lat))
    latmax = np.ceil(np.max(lat))
    plt.figure(figsize=(10,7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())
    if type(bgimage) != str:
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90], 
                   zorder=1)
    #end if
    vmin = 0
    vmax = 1
    im2 = plt.pcolormesh(lon, lat, values[:,:,0], 
                         cmap='hot_r', vmin=vmin, 
                         vmax=vmax, zorder=2)
    ax.coastlines(zorder=4)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    plt.colorbar(im2, ax=ax, label='Normed probability density')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('2D normed spatial PDF map')
    name = obj.stat_obj.name
    path = os.path.join(wd_images, name)
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    filename = os.path.join(path, '2D_spatial_PDF_map.png')
    plt.savefig(filename)
    plt.show()
    return
#end func
    
def plot_coverage_region(domain_dict, path, bgimage=""):
    """
    Plots the region of coverage of a geophysical model, both with and without
    areas where deposits are able to appear at the surface.
    
    
    Parameter
    ---------
    domain_dict : dictionary
        Dictionary of polygons which must have at least 'region' and
        'uncovered_region' polygons within it.
        
    path : string
        Directory in which to save image.
        
    bgimage : string or matplotlib.imread object
        Background image on which to overlay plot.
    
    
    """
    region = domain_dict['region']
    allowed_region = domain_dict['allowed_region']
    minlon, minlat, maxlon, maxlat = \
        np.array(region.bounds) + np.array([-0.5, -0.5, 0.5, 0.5])
    minlon = np.max([minlon, -180])
    maxlon = np.min([maxlon, 180])
    minlat = np.max([minlat, -90])
    maxlat = np.min([maxlat, 90])
    
    bounds = np.array([minlon, maxlon, minlat, maxlat])
    
    plt.figure(figsize=(10,7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())
    if type(bgimage) != str:
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90])
    #end if
    ax.coastlines(zorder=5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent(bounds)
    if type(region) == Polygon:
        xs, ys = region.exterior.xy    
        ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
    else:
        for geom in region.geoms: 
            if type(geom) == Point or type(geom) == LineString:
                continue
            #end if
            xs, ys = geom.exterior.xy    
            ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
        #end for
    #end if
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Coverage')
    plt.savefig(os.path.join(path, 'domain.png'))
    plt.show()
    
    plt.figure(figsize=(10,7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())
    if type(bgimage) != str:
        plt.imshow(bgimage, 
                   origin='upper', 
                   transform=ccrs.PlateCarree(), 
                   extent=[-180, 180, -90, 90])
    #end if
    ax.coastlines(zorder=5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent(bounds)
    if type(region) == Polygon:
        xs, ys = region.exterior.xy    
        ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
    else:
        for geom in region.geoms:    
            if type(geom) == Point or type(geom) == LineString:
                continue
            #end if
            xs, ys = geom.exterior.xy    
            ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
        #end for
    #end if

    if type(allowed_region) == Polygon:
        xs, ys = allowed_region.exterior.xy
        ax.fill(xs, ys, alpha=0.5, fc='b', ec='none')
    else:
        for geom in allowed_region.geoms:
            if type(geom) == Point or type(geom) == LineString:
                continue
            #end if
            xs, ys = geom.exterior.xy
            ax.fill(xs, ys, alpha=0.5, fc='b', ec='none')
        #end for
    # end if
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Coverage')
    plt.savefig(os.path.join(path, 'uncovered_domain.png'))
    plt.show()
#end func