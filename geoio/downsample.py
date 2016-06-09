import numpy as np
import itertools
import logging
from numba import jit, vectorize
from math import ceil, floor

logger = logging.getLogger(__name__)

@jit(nopython=True)
def test1():
    a=1
    return a+1

@vectorize
def test2():
    a = np.array([1,2,3])
    return a.sum()


# # Need to debug the following:
# Out[204]:
# array([[[ 161.25,  331.75,  389.75,  328.5 ,  372.75],
#         [ 204.75,  432.25,  417.25,  336.  ,  348.  ],
#         [ 215.  ,  424.5 ,  398.25,  311.75,  326.25],
#         [ 202.  ,  348.  ,  309.  ,  300.5 ,  340.5 ],
#         [ 204.5 ,  319.75,  307.75,  310.75,  311.  ]]])
#
# In [205]: geoio.utils.block_view(data[0,:,:].squeeze(),block=(2,2),strides=(2,2)).mean(axis=(2,3))
# Out[205]:
# array([[ 161.25,  204.75,  215.  ,  202.  ,  204.5 ],
#        [ 331.75,  432.25,  424.5 ,  348.  ,  319.75],
#        [ 389.75,  417.25,  398.25,  309.  ,  307.75],
#        [ 328.5 ,  336.  ,  311.75,  300.5 ,  310.75],
#        [ 372.75,  348.  ,  326.25,  340.5 ,  311.  ]])
#
# #### These are opposite!!!! ####



def downsample(arr,
               shape = None,
               factor = None,
               ul_corner = None,
               lr_corner = None,
               method = 'aggregate'):

    # If arr comes in as a 2D array, assume this is a single band from
    # a 3D image and reshape accordingly
    if len(arr.shape) == 2:
        arr = arr[np.newaxis, :, :]

    if shape and factor:
        raise ValueError('Either shape or factor can be specificed, not both.')

    if not shape and not factor:
        raise ValueError('Either shape or factor needs to be specified.')

    if method not in ['aggregate','nearest']:
        raise ValueError("The downsample method can be 'aggregate' or "
                         "'nearest'.")

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

    # Set x_steps and y_steps for the downsampling process below
    if shape:
        x_start = 0
        x_stop = arr.shape[1] # no -1 to bracket for base 0 indexing
        x_num = shape[0]+1 # +1 to index on both sides of block.
        y_start = 0
        y_stop = arr.shape[2] # no -1 to bracket for base 0 indexing
        y_num = shape[1]+1 # +1 to index on both sides of block.

    if factor:
        x_start = 0
        x_stop = arr.shape[1]
        x_num = int(round(arr.shape[1]*factor[0]))+1
        y_start = 0
        y_stop = arr.shape[2]
        y_num = int(round(arr.shape[2]*factor[1]))+1

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

    print(x_steps)
    print(y_steps)

    # import ipdb; ipdb.set_trace()

    if method == 'aggregate':
        return aggregate(arr,x_steps,y_steps)
    elif method == 'nearest':
        raise NotImplementedError()

#@jit(nopython=True)
def aggregate(arr,x_steps,y_steps):

    # x_index = range(len(x_steps)-1)
    # y_index = range(len(y_steps)-1)

    out = np.empty((arr.shape[0],len(x_steps)-1, len(y_steps)-1))


    # for b in xrange(out.shape[0]):
    #     for x,y in itertools.product(x_index,y_index):
    #         # inner
    #         x_inner = np.arange(np.ceil(x_steps[x]),np.floor(x_steps[x+1])+1)
    #         y_inner = np.arange(np.ceil(y_steps[y]),np.floor(y_steps[y+1])+1)
    #         agg = 0
    #         out[b,x,y] = arr[b,x_index[x]:x_index[x+1],y_index[y]:y_index[y+1]]

    if logger.isEnabledFor(logging.DEBUG):
        talk = True
    else:
        talk = False

    for b in xrange(out.shape[0]):
        for x in xrange(out.shape[1]):
            for y in xrange(out.shape[2]):

                print(x)
                print(y)

                # import ipdb; ipdb.set_trace()

                # Debug info for initial state of this block
                if talk:
                    logger.debug('*** Starting new block calculation ***')
                    logger.debug('Block steps left/right: %s, %s' %
                                                        (x_steps[x], x_steps[x+1]))
                    logger.debug('Block steps top/bottom: %s, %s' %
                                                        (y_steps[y], y_steps[y+1]))

                # initialize sum variable
                s = 0

                # sum center pixels
                left = int(ceil(x_steps[x]))
                right = int(floor(x_steps[x+1]))
                top = int(ceil(y_steps[y]))
                bottom =  int(floor(y_steps[y+1]))
                # if talk:
                logger.debug('left, right, top, bottom: %s, %s, %s, %s' %
                                                (left, right, top, bottom))

                # import ipdb; ipdb.set_trace()

                # for p in xrange(left,right+1,1):
                #     for q in xrange(bottom,top+1,1):
                #         logger.debug(p,q)
                #         logger.debug(arr[b,p,q])
                #         s += arr[b,p,q]
                #         logger.debug(s)
                #         logger.debug('---')

                s += arr[b,left:right,top:bottom].sum()

                # sum edges

                # sum corners

                # calculate weight
                weight = (x_steps[x+1]-x_steps[x])*(y_steps[y+1]-y_steps[y])

                if talk:
                    logger.debug('sum is: %s' % s)
                    logger.debug('weight is: %s' % weight)

                    logger.debug('*** block sum complete ***')

                out[b,y,x] = s/float(weight)

    return out
