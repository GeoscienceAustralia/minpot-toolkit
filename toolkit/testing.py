"""
Boilerplate
"""

import numpy as np
import os
from toolkit.stat_corr import compute_normed_cdf
from toolkit.main import MineralPotential #Used to unpickle save files.

def check_cdf_equivalence(cdf_array1, cdf_array2):
    """
    Checks to see if the cdf array produced by the updated toolkit remains the
    same as that produced by the original toolkit.
    
    
    """
    l = len(cdf_array1)
    l2 = len(cdf_array2)
    if l != l2: return False
    for i in range(l):
        m = max(1.0, max(cdf_array1[i]['cdf']))
        m2 = max(1.0, max(cdf_array2[i]['cdf']))
        if m != m2: return False
        diff = cdf_array1[i]['cdf'] - cdf_array2[i]['cdf']
        chk = (diff > 0.05*m)
        if np.any(chk): return False
        cdfr1 = np.mean(cdf_array1['cdf_random'][i],axis=0)
        cdfr2 = np.mean(cdf_array2['cdf_random'][i],axis=0)
        cdfstd1 = np.std(cdf_array1['cdf_random'][i],axis=0)
        cdfstd2 = np.std(cdf_array2['cdf_random'][i],axis=0)
        diff = np.abs(cdfr1 - cdfr2)
        std = np.max([cdfstd1, cdfstd2], axis=0)
        chk = (np.abs(diff)>std)
        if np.any(chk): return False
    #end for
    return True
#end func

def check_D_value_equivalence(dvalues1, dvalues2):
    """
    Checks to see if the D values produced by the updated toolkit remains the
    same as that produced by the original toolkit.
    
    
    """
    l = len(dvalues1)
    l2 = len(dvalues2)
    if l != l2: return False
    diff = dvalues1 - dvalues2
    avg = (dvalues1 + dvalues2)/2.0
    ind = (avg != 0.0)
    chk = np.all(np.abs(diff[ind])/avg[ind] < 0.05)
    return chk
#end func

def check_distance_array_equivalence(distance_array1, distance_array2):
    """
    Checks to see if the distance array produced by the updated toolkit remains 
    the same as that produced by the original toolkit.
    
    
    """
    l = len(distance_array1)
    l2 = len(distance_array2)
    if l != l2: return False
    for i in range(l):
        #check which coordinates of array NaN is present
        chk1 = np.logical_and(np.isnan(distance_array1[i]['distance']),
                              np.isnan(distance_array2[i]['distance']))
        #check if non NaN values are equivalent
        diff = distance_array1[i]['distance'] - distance_array2[i]['distance']
        avg = (distance_array1[i]['distance'] + \
               distance_array2[i]['distance'])/2
        chk2 = (diff/avg < 0.01)
        #return truth value at all coordinates
        chk3 = np.logical_xor(chk1, chk2)
        #check if all true
        if ~np.all(chk3): return False
    #end for
    return True
#end func

def check_P_value_equivalence(kolmogorov_smirnov1, kolmogorov_smirnov2):
    """
    Checks to see if the Kolmogorov-Smirnov test values produced by the updated 
    toolkit remains the same as that produced by the original toolkit.
    
    
    """
    l = len(kolmogorov_smirnov1)
    l2 = len(kolmogorov_smirnov2)
    if l != l2: return False
    diff = kolmogorov_smirnov1 - kolmogorov_smirnov2
    avg = (kolmogorov_smirnov1 + kolmogorov_smirnov2)/2.0
    chk = np.all(np.abs(diff)/avg < 0.05)
    return chk
#end func
    
def process(path1, path2, testname):
    """
    Perform tests to check if results remain consistend using new toolkit over
    old toolkit.
    
    
    Calls
    -----
    check_p_value_equivalence
    
    check_D_value_equivalence
    
    check_cdf_array_equivalence
    
    check_distance_array_equivalence
    
    
    """
    savefile = os.path.join(path1, 'savefile.npy')
    cdf2file = os.path.join(path2, 'cdf.npy')
    dist2file = os.path.join(path2, 'test_points_100ohmm_distance_array.npy')
    
    obj = np.load(savefile, allow_pickle=True).item()
    cdf_obj1 = obj.StatsManager.stat_obj_dict[testname]
    n_deposits = cdf_obj1.n_deposits
    
    dist_array2 = np.load(dist2file, allow_pickle=True)
    cdf_array2 = np.load(cdf2file, allow_pickle=True)
    cdf_normed2, cdf_random_normed2, cdf_random_std2, dvalues2, dmax2 = \
        compute_normed_cdf(cdf_array2, n_deposits, 'Positive')
    kolmogorov_smirnov2 = np.exp(-2.*(n_deposits**2)*(dmax2**2)/(2*n_deposits))
    
    t1 = check_cdf_equivalence(cdf_obj1.cdf_array, cdf_array2)
    t2 = check_D_value_equivalence(cdf_obj1.dvalues, dvalues2)
    t3 = check_distance_array_equivalence(cdf_obj1.distance_array, dist_array2)
    t4 = check_P_value_equivalence(cdf_obj1.kolmogorov_smirnov, 
                                   kolmogorov_smirnov2)
    return t1, t2, t3, t4
    
if __name__ == '__main__':
    path1 = r'C:/Users/U37509/Documents/Mineral_Potential/outputs/'
    path2 = r'C:/Users/U37509/Documents/Mineral_Potential/mineral-potential-toolkit-original/data/outputs/'
    testname = 'test17'
    tests = process(path1, path2, testname)
#end if
