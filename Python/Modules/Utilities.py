# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:46:30 2016

@author: Alonyan
"""
import numpy as np
import scipy.io as sio

def num2str(num, precision): 
    return "%0.*f" % (precision, num)
    
def colorcode(datax, datay):
    from scipy import interpolate
    import numpy as np
    H, xedges, yedges = np.histogram2d(datax,datay, bins=30)
    xedges = (xedges[:-1]+xedges[1:])/2
    yedges = (yedges[:-1]+yedges[1:])/2
    f = interpolate.RectBivariateSpline(xedges,yedges , H)
    
    z = np.array([])
    for i in datax.index:
        z = np.append(z,f(datax[i],datay[i]))
    #z=(z-min(z))/(max(z)-min(z))   
    z[z<0] = 0
    idx = z.argsort()
    return z, idx
    
class kmeans:        
    def __init__(self, X, K):
        # Initialize to K random centers
        oldmu = X.sample(K).values#np.random.sample(X, K)
        mu = X.sample(K).values#np.random.sample(X, K)
        while not _has_converged(mu, oldmu):
            oldmu = mu
            # Assign all points in X to clusters
            clusters = _cluster_points(X, mu)
            # Reevaluate centers
            mu = _reevaluate_centers(oldmu, clusters)
        self.mu = mu
        self.clusters = clusters
        #return(mu, clusters)
        
def _cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def _reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def _has_converged(mu, oldmu):
    return (set(mu) == set(oldmu))
    
    
    
    
def makeTicks():
    a = np.outer(np.arange(1,10),10**np.arange(1,2)).T.reshape((1,-1)).squeeze()
    ticks = np.append(-a[::-1],0)
    ticks = np.append(-100,ticks)
    a = np.outer(np.arange(1,10),10**np.arange(1,6)).T.reshape((1,-1)).squeeze()
    ticks = np.append(ticks,a[:])
    emptvec = ['','','','','','','','']
    ticklabels = ['-0.1']+emptvec+['']+['0']+emptvec+['']+['0.1']+emptvec+['1']+emptvec+['10']+emptvec+['100']+emptvec
    return ticks, ticklabels    
    
    
    

#Utils for opening MAT files
def print_mat_nested(d, indent=0, nkeys=0):
    """Pretty print nested structures from .mat files   
    Inspired by: `StackOverflow <http://stackoverflow.com/questions/3229419/pretty-printing-nested-dictionaries-in-python>`_
    """
    
    # Subset dictionary to limit keys to print.  Only works on first level
    if nkeys>0:
        d = {k: d[k] for k in d.keys()[:nkeys]}  # Dictionary comprehension: limit to first nkeys keys.

    if isinstance(d, dict):
        for key, value in d.iteritems():         # iteritems loops through key, value pairs
          print '\t' * indent + 'Key: ' + str(key)
          print_mat_nested(value, indent+1)

    if isinstance(d,np.ndarray) and d.dtype.names is not None:  # Note: and short-circuits by default
        for n in d.dtype.names:    # This means it's a struct, it's bit of a kludge test.
            print '\t' * indent + 'Field: ' + str(n)
            print_mat_nested(d[n], indent+1)


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    
    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    dict1 = {}
    for key in dict:
        if isinstance(dict[key],np.ndarray):
            i=1
            for inst in dict[key]:
                if isinstance(inst, sio.matlab.mio5_params.mat_struct):
                    dict1[key+'_'+str(i)] = _todict(inst) 
                    i+=1
        elif isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict1[key] = _todict(dict[key])
            
    return dict1       

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        elif isinstance(elem,np.ndarray):
            dict[strg] = _tolist(elem)
        else:
            dict[strg] = elem
    return dict

def _tolist(ndarray):
    '''
    A recursive function which constructs lists from cellarrays 
    (which are loaded as numpy ndarrays), recursing into the elements
    if they contain matobjects.
    '''
    elem_list = []            
    for sub_elem in ndarray:
        if isinstance(sub_elem, sio.matlab.mio5_params.mat_struct):
            elem_list.append(_todict(sub_elem))
        elif isinstance(sub_elem,np.ndarray):
            elem_list.append(_tolist(sub_elem))
        else:
            elem_list.append(sub_elem)
    return elem_list
#
#def _todict(matobj):
#    '''
#    A recursive function which constructs from matobjects nested dictionaries
#    '''
#    dict = {}
#    for strg in matobj._fieldnames:
#        elem = matobj.__dict__[strg]
#        if isinstance(elem, np.ndarray):
#            i=1
#            for el in elem:
#                if isinstance(el, sio.matlab.mio5_params.mat_struct):
#                    dict[strg+'_'+str(i)] = _todict(el) 
#                    i+=1
#                else:
#                    dict[strg] = elem
#        elif isinstance(elem, sio.matlab.mio5_params.mat_struct):
#            dict[strg] = _todict(elem)
#        else:
#            dict[strg] = elem
#    return dict
#
#    