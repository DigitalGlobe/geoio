import numpy as np
from math import ceil, floor
from numba import (jit, guvectorize, vectorize, float64, float32,
                   int16, uint16, double)

__all__ = ['aggregate_pixel', 'nearest_pixel', 'max_pixel', 'min_pixel', 'aggregate_numba_3d', 'aggregate_guvec', 'nearest_guvec',
           'max_guvec', 'min_guvec']

@jit(nopython=True)
def aggregate_pixel(arr,x_step,y_step):
    """Aggregation code for a single pixel"""

    # Set x/y to zero to mimic the setting in a loop
    # Assumes x_step and y_step in an array-type of length 2
    x = 0
    y = 0

    # initialize sum variable
    s = 0.0

    # sum center pixels
    left = int(ceil(x_step[x]))
    right = int(floor(x_step[x+1]))
    top = int(ceil(y_step[y]))
    bottom =  int(floor(y_step[y+1]))
    s += arr[left:right,top:bottom].sum()

    # Find edge weights
    wl = left - x_step[x]
    wr = x_step[x+1] - right
    wt = top - y_step[y]
    wb = y_step[y+1] - bottom
    # sum edges - left
    s += arr[left-1:left,top:bottom].sum() * wl
    # sum edges - right
    s += arr[right:right+1,top:bottom].sum() * wr
    # sum edges - top
    s += arr[left:right,top-1:top].sum() * wt
    # sum edges - bottom
    s += arr[left:right,bottom:bottom+1].sum() * wb

    # sum corners ...
    # ul
    s += arr[left-1:left,top-1:top].sum() * wl * wt
    # ur
    s += arr[right:right+1,top-1:top].sum() * wr * wt
    # ll
    s += arr[left-1:left,bottom:bottom+1].sum() * wl * wb
    # lr
    s += arr[right:right+1,bottom:bottom+1].sum() * wr * wb

    # calculate weight
    weight = (x_step[x+1]-x_step[x])*(y_step[y+1]-y_step[y])

    return s/float(weight)

@jit(nopython=True)
def nearest_pixel(arr,x_step,y_step):
    """Aggregation code for a single pixel"""

    # Set x/y to zero to mimic the setting in a loop
    # Assumes x_step and y_step in an array-type of length 2
    x = 0
    y = 0

    # initialize sum variable
    s = 0.0

    # nearest neighbor
    x_center = int(np.mean(x_step[x:x+2]))
    y_center = int(np.mean(x_step[y:y+2]))
    s += arr[x_center,y_center]

    return s

@jit(nopython=True)
def max_pixel(arr,x_step,y_step):
    """Aggregation code for a single pixel"""

    # Set x/y to zero to mimic the setting in a loop
    # Assumes x_step and y_step in an array-type of length 2
    x = 0
    y = 0

    # initialize sum variable
    s = 0.0

    # sum center pixels
    left = int(ceil(x_step[x]))
    right = int(floor(x_step[x+1]))
    top = int(ceil(y_step[y]))
    bottom =  int(floor(y_step[y+1]))
    s += arr[left-1:right+1,top-1:bottom+1].max()

    return s

@jit(nopython=True)
def min_pixel(arr,x_step,y_step):
    """Aggregation code for a single pixel"""

    # Set x/y to zero to mimic the setting in a loop
    # Assumes x_step and y_step in an array-type of length 2
    x = 0
    y = 0

    # initialize sum variable
    s = 0.0

    # sum center pixels
    left = int(ceil(x_step[x]))
    right = int(floor(x_step[x+1]))
    top = int(ceil(y_step[y]))
    bottom =  int(floor(y_step[y+1]))
    s += arr[left-1:right+1,top-1:bottom+1].min()

    return s


@jit(nopython=True)
def aggregate_numba_3d(arr,x_steps,y_steps):
    """TBD"""
    out = np.zeros((arr.shape[0],len(x_steps)-1,len(y_steps)-1),
                   dtype=arr.dtype)

    for b in xrange(out.shape[0]):
        for x in xrange(out.shape[1]):
            for y in xrange(out.shape[2]):
                out[b,x,y] = aggregate_pixel(arr[b,:,:],
                                             x_steps[x:x+2],
                                             y_steps[y:y+2])

    return out

# The types handled are the same in contstants.py DICT_GDAL_TO_NP
@guvectorize(['void(uint8[:,:],float64[:],float64[:],uint8[:,:])',
              'void(uint16[:,:],float64[:],float64[:],uint16[:,:])',
              'void(uint32[:,:],float64[:],float64[:],uint32[:,:])',
              'void(int16[:,:],float64[:],float64[:],int16[:,:])',
              'void(int32[:,:],float64[:],float64[:],int32[:,:])',
              'void(float32[:,:],float64[:],float64[:],float32[:,:])',
              'void(float64[:,:],float64[:],float64[:],float64[:,:])'],
            '(a,b),(c),(d),(m,n)',target='parallel',nopython=True)
def aggregate_guvec(arr, x_steps, y_steps, out):
    """TBD"""
    for x in xrange(out.shape[0]):
        for y in xrange(out.shape[1]):
            out[x,y] = aggregate_pixel(arr,x_steps[x:x+2],y_steps[y:y+2])

# The types handled are the same in contstants.py DICT_GDAL_TO_NP
@guvectorize(['void(uint8[:,:],float64[:],float64[:],uint8[:,:])',
              'void(uint16[:,:],float64[:],float64[:],uint16[:,:])',
              'void(uint32[:,:],float64[:],float64[:],uint32[:,:])',
              'void(int16[:,:],float64[:],float64[:],int16[:,:])',
              'void(int32[:,:],float64[:],float64[:],int32[:,:])',
              'void(float32[:,:],float64[:],float64[:],float32[:,:])',
              'void(float64[:,:],float64[:],float64[:],float64[:,:])'],
             '(a,b),(c),(d),(m,n)', target='parallel',
             nopython=True)
def nearest_guvec(arr, x_steps, y_steps, out):
    """TBD"""
    for x in xrange(out.shape[0]):
        for y in xrange(out.shape[1]):
            out[x, y] = nearest_pixel(arr,x_steps[x:x + 2],y_steps[y:y + 2])

# The types handled are the same in contstants.py DICT_GDAL_TO_NP
@guvectorize(['void(uint8[:,:],float64[:],float64[:],uint8[:,:])',
              'void(uint16[:,:],float64[:],float64[:],uint16[:,:])',
              'void(uint32[:,:],float64[:],float64[:],uint32[:,:])',
              'void(int16[:,:],float64[:],float64[:],int16[:,:])',
              'void(int32[:,:],float64[:],float64[:],int32[:,:])',
              'void(float32[:,:],float64[:],float64[:],float32[:,:])',
              'void(float64[:,:],float64[:],float64[:],float64[:,:])'],
             '(a,b),(c),(d),(m,n)', target='parallel',
             nopython=True)
def max_guvec(arr, x_steps, y_steps, out):
    """TBD"""
    for x in xrange(out.shape[0]):
        for y in xrange(out.shape[1]):
            out[x, y] = max_pixel(arr,x_steps[x:x + 2],y_steps[y:y + 2])


# The types handled are the same in contstants.py DICT_GDAL_TO_NP
@guvectorize(['void(uint8[:,:],float64[:],float64[:],uint8[:,:])',
              'void(uint16[:,:],float64[:],float64[:],uint16[:,:])',
              'void(uint32[:,:],float64[:],float64[:],uint32[:,:])',
              'void(int16[:,:],float64[:],float64[:],int16[:,:])',
              'void(int32[:,:],float64[:],float64[:],int32[:,:])',
              'void(float32[:,:],float64[:],float64[:],float32[:,:])',
              'void(float64[:,:],float64[:],float64[:],float64[:,:])'],
             '(a,b),(c),(d),(m,n)', target='parallel',
             nopython=True)
def min_guvec(arr, x_steps, y_steps, out):
    """TBD"""
    for x in xrange(out.shape[0]):
        for y in xrange(out.shape[1]):
            out[x, y] = min_pixel(arr,x_steps[x:x + 2],y_steps[y:y + 2])