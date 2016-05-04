import tinytools as tt
import numpy as np
import matplotlib.pyplot as plt

# Currently doesn't return figure handles to help with memory size issues on
# large images.  This could probably use a close look at ways to make it
# more efficient:
# https://pypi.python.org/pypi/ModestImage
# https://github.com/ipython/ipython/issues/1623/
# http://stackoverflow.com/questions/15345336/memory-leak-in-matplotlib-imshow

def imshow(data,stretch=[0.02,0.98],stretch_type='linear'):
    """Convenience method to do all the plotting gymnastics to get a resonable
    looking image plot.

    Input:
    data            numpy array in gdal band order - 3 dimensions (bands, x, y)
    stretch         stretch values on a scale of [0,1]
    stretch_type    type of stretch scale (only linear is curretly supported)
    """

    if len(data[:,0,0]) != 3:
        raise ValueError('This convenience function is only implemented ' \
                         'for three bands.  Use img.get_data(bands=...) to ' \
                         'retrieve specific data.')
                        # ToDo - This could also be speed up by indexing a
                        # single numpy array.

    # define stretch
    # Possibly useful code for additional stretches at:
    # http://scikit-image.org/docs/dev/api/skimage.exposure.html
    # also
    # http://scikit-image.org/docs/dev/auto_examples/plot_equalize.html
    if stretch_type == "linear":
        pass
    else:
        raise ValueError('The passed value of stretch is not implemented.')

    # Get the per-band scaled data
    data = tt.np_img.conv_to_bandslast(data)
    data = data.astype('float32')
    lims = np.percentile(data,(2,98),axis=(0,1))
    for x in xrange(len(data[0,0,:])):
        top = lims[:,x][1]
        bottom = lims[:,x][0]
        data[:,:,x] = (data[:,:,x]-bottom)/float(top-bottom)
    data = np.clip(data,0,1)

    plt.imshow(data);
    plt.show(block=False)

    # ToDo: fix the handle return for later update and memory issues - see
    # comments above
    # return handle

def hist(data):
    """Convenience method to do all the plotting gymnastics to get a resonable
    looking histogram plot.

    Input: data - numpy array in gdal format - (bands, x, y)

    Returns:  matplotlib figure handle

    Adapted from:
    http://nbviewer.jupyter.org/github/HyperionAnalytics/PyDataNYC2014/blob/master/color_image_processing.ipynb
    """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.hold(True)
    for x in xrange(len(data[:,0,0])):
        counts, edges = np.histogram(data[x,:,:],bins=100)
        centers = [(edges[i]+edges[i+1])/2.0 for i,v in enumerate(edges[:-1])]
        ax1.plot(centers,counts)
    plt.hold(False)

    plt.show(block=False)

    # return fig
