import logging
import numpy as np

### Optional imports ###
try:
    import cv2
    use_cv2 = True
except ImportError:
    use_cv2 = False

try:
    from numba import (jit, guvectorize, vectorize, float64, float32,
                       int16, uint16, double)
    from downsample_numba import *
    use_numba = True
except ImportError:
    use_numba = False

logger = logging.getLogger(__name__)

def downsample(arr,
               shape = None,
               factor = None,
               extent = None,
               method = 'aggregate',
               no_data_value = None,
               source = None):
    """
    Downsampling an array with either a custom routine from below or
    cv2 based on which options is available for a given downsample.

    Parameters
    ----------
    arr : array_like
        Image data in a two or three dimension array in band-first order.
    shape : list
        Shape of the desired output array.
    factor : integer, float, or length two iterable
        Factor by which to scale the image (must be less than one).
    extent : length four array-like
        List of upper-left and lower-right corner coordinates for the edges
        of the resample.  Coordinates should be expressed in pixel space.
        i.e. [-1,-2,502,501]
    method : strings
        Method to use for the downsample - 'aggregate', 'nearest', 'max', or
        'min.
    no_data_value : int
        Data value to treat as no_data
    source : strings
        Package to use for algorithm - opencv ('cv2') or custom 'numba' below.

    Returns
    -------
    ndarray
        Three dimensional numpy array of downsampled data
    """

    # If arr comes in as a 2D array, assume this is a single band from
    # a 3D image and reshape accordingly
    if len(arr.shape) == 2:
        arr = arr[np.newaxis, :, :]

    # Check input parameters
    if shape and factor:
        raise ValueError('Either shape or factor can be specificed, not both.')

    if not shape and not factor:
        raise ValueError('Either shape or factor needs to be specified.')

    if extent is not None:
        if (len(extent) != 4):
            raise ValueError('extent needs to be an array-like object of '
                             'length four.  It should desribe the upper-left '
                             'and lower-right corners of '
                             'the desired resmple area in pixels space as '
                             '[ul_x, ul_y, lr_x, lr_y] and can describe '
                             'points outside the image extent, '
                             'i.e. [-1, -2.1, 500, 501].  The shape or '
                             'factor parameter must also be provided to '
                             'describe the size of the requested array.')

    if factor is not None:
        # Prep factor based on input type/format
        if isinstance(factor,(float,int)):
            factor = [factor,factor]

        # Check that factor values are valid
        for f in factor:
            if f >= 1:
                raise ValueError('Factor values should be less than one.')

    if shape is not None:
        if (shape[0] >= arr.shape[1]) | (shape[1] >= arr.shape[2]):
            raise ValueError('The requested downsample shape should be less '
                             'than the array passed in.')

    # Check other input parameters
    if method not in ['aggregate','nearest','max','min']:
        raise ValueError("The downsample method can be 'aggregate' or "
                         "'nearest'.")

    # Set x_steps and y_steps for the downsampling process below
    if shape is not None:
        x_start = 0
        x_stop = arr.shape[1] # no -1 to bracket for base 0 indexing
        x_num = shape[0]+1 # +1 to index on both sides of block.
        y_start = 0
        y_stop = arr.shape[2] # no -1 to bracket for base 0 indexing
        y_num = shape[1]+1 # +1 to index on both sides of block.

    if factor is not None:
        x_start = 0
        x_stop = arr.shape[1]
        x_num = int(round(arr.shape[1]*factor[0]))+1
        y_start = 0
        y_stop = arr.shape[2]
        y_num = int(round(arr.shape[2]*factor[1]))+1

    if extent is not None:
        x_start = extent[0]
        x_stop = extent[2]
        y_start = extent[1]
        y_stop = extent[3]

    x_steps = np.linspace(x_start,x_stop,x_num)
    y_steps = np.linspace(y_start,y_stop,y_num)

    logger.debug('x_start, x_stop, x_num: %s, %s, %s' %
                                                (x_start, x_stop, x_num))
    logger.debug('y_start, y_stop, y_num: %s, %s, %s' %
                                                (y_start, y_stop, y_num))
    logger.debug('length of x_steps: %s' % len(x_steps))
    logger.debug('length of y_steps: %s' % len(y_steps))

    logger.debug('beggining of x_steps: %s ...' % x_steps[:3])
    logger.debug('end of x_steps: ... %s' % x_steps[-3:])
    logger.debug('beggining of y_steps: %s ...' % y_steps[:3])
    logger.debug('end of y_steps: ... %s' % y_steps[-3:])

    return downsample_to_grid(arr,x_steps,y_steps,no_data_value,method,source)

def downsample_to_grid(arr,x_steps,y_steps,no_data_value=None,
                                            method='aggregate',source=None):
    """
    Function to choose which execuatable to use based on efficiency and
    availability.

    Parameters
    ----------
    arr
    x_steps
    y_steps
    method

    Returns
    -------

    """

    global use_cv2
    global use_numba
    if source is None:
        pass
    elif source=='cv2':
        use_cv2 = True
        use_numba = False
    elif source=='numba':
        use_cv2 = False
        use_numba = True
    else:
        raise ValueError('Invalid source provided.')

    if no_data_value:
        arr = np.where(arr == no_data_value, 0, arr)

    # Find the appropriate algorithm to run with.
    if ((x_steps[0] == 0) and (x_steps[-1] == arr.shape[1]) and
        (y_steps[0] == 0) and (y_steps[-1] == arr.shape[2]) and
        use_cv2):
        # If the requested steps are exact subsets of the image array and
        # cv2 is available, use cv2 since it is fast.
        if method == 'aggregate':
            logger.debug('running aggregate with opencv::resize')
            type_cv_code = cv2.INTER_AREA
            out = run_opencv_resize(arr,x_steps,y_steps,type_cv_code)
        elif method == 'nearest':
            logger.debug('running nearest neighbor downsample with '
                         'opencv::resize')
            type_cv_code = cv2.INTER_NEAREST
            out = run_opencv_resize(arr,x_steps,y_steps,type_cv_code)
        else:
            raise ValueError('Specified method is not available')
    elif use_numba:
        # If cv2 isn't available or the requested steps are from a grid
        # that doesn't nicely overlap this image, use custom implementations
        # from below.
        if method == 'aggregate':
            logger.debug('running aggregate with custom numba function.')
            out = run_numba_aggregate(arr,x_steps,y_steps)
        elif method == 'nearest':
            logger.debug('running nearest neighbor downsample with '
                         'custom numba function.')
            out = run_numba_nearest(arr, x_steps, y_steps)
        elif method == 'max':
            logger.debug('running max with custom numba function.')
            out = run_numba_max(arr, x_steps, y_steps)
        elif method == 'min':
            logger.debug('running min with custom numba function.')
            out = run_numba_min(arr, x_steps, y_steps)
        else:
            raise ValueError('Specified method is not available')
    else:
        raise ValueError('No downsampling routine available for the '
                         'requested parameters.  Either opencv or numba are '
                         'needed to run the downsample calculations.  '
                         'Additionally, if the requested grid does not '
                         'align with the edges of the image, you can not '
                         'use opencv and will need numba to run this '
                         'function.  You can always just use get_data() '
                         'and use an external resampling routine!')

    if no_data_value is not None:
        return np.where(out == 0, no_data_value, out)
    else:
        return out

def run_opencv_resize(arr,x_steps,y_steps,type_cv_code):
    """TBD"""
    size = (len(x_steps)-1,len(y_steps)-1)
    out = np.empty([arr.shape[0]]+list(size))

    for b in xrange(out.shape[0]):
        out[b,:,:] = cv2.resize(arr[b,:,:],dsize=size[::-1],dst=None,
                                                interpolation=type_cv_code)
    return out

def run_numba_aggregate(arr,x_steps,y_steps):
    """TBD"""
    use_jit = False

    if use_jit:
        out = aggregate_numba_3d(arr,x_steps,y_steps)
    else:
        out = np.zeros((arr.shape[0],len(x_steps)-1,len(y_steps)-1),
                       dtype=arr.dtype)
        aggregate_guvec(arr,x_steps,y_steps,out)

    return out

def run_numba_nearest(arr,x_steps,y_steps):
    """TBD"""
    out = np.zeros((arr.shape[0],len(x_steps)-1,len(y_steps)-1),
                   dtype=arr.dtype)
    nearest_guvec(arr,x_steps,y_steps,out)

    return out

def run_numba_max(arr,x_steps,y_steps):
    """TBD"""
    out = np.zeros((arr.shape[0],len(x_steps)-1,len(y_steps)-1),
                   dtype=arr.dtype)
    max_guvec(arr,x_steps,y_steps,out)

    return out

def run_numba_min(arr,x_steps,y_steps):
    """TBD"""
    out = np.zeros((arr.shape[0],len(x_steps)-1,len(y_steps)-1),
                   dtype=arr.dtype)
    min_guvec(arr,x_steps,y_steps,out)

    return out


def main():
    import time
    import dgsamples
    import geoio

    img_small = geoio.GeoImage(dgsamples.wv2_longmont_1k.ms)
    data_small = img_small.get_data()

    start = time.time()
    out_small_numba = downsample(data_small,shape=[300,300],source='numba')
    print('small numba:  %s' % (time.time()-start))

    start = time.time()
    out_small_cv2 = downsample(arr=data_small,shape=(300,300),source='cv2')
    print('small cv2:  %s' % (time.time()-start))

    print('Max diff is:  %s' % (out_small_numba-out_small_cv2).max())
    print('Min diff is:  %s' % (out_small_numba-out_small_cv2).min())

    # img_big = geoio.GeoImage('/mnt/panasas/nwl/data_HIRES/Gibraltar/VNIR/054312817010_01_P001_MUL/15FEB28112650-M2AS_R1C1-054312817010_01_P001.TIF')
    # data_big = img_big.get_data()
    #
    # start = time.time()
    # out_big_numba = downsample(data_big,shape=[1000,1000],source='numba')
    # print('big numba:  %s' % (time.time()-start))
    #
    # start = time.time()
    # out_big_cv2 = downsample(arr=data_big,shape=(1000,1000),source='cv2')
    # print('big cv2:  %s' % (time.time()-start))
    #
    # print('Max diff is:  %s' % (out_big_numba-out_big_cv2).max())
    # print('Min diff is:  %s' % (out_big_numba-out_big_cv2).min())

if __name__ == "__main__":
    main()