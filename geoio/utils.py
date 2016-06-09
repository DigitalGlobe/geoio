from numpy.lib.stride_tricks import as_strided

def block_view(A, block=(3,3), strides=(2,2)):
    """Provide a 2D block view to 2D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved ;-)
    shape= (A.shape[0]/ block[0], A.shape[1]/ block[1])+ block
    strides= (strides[0]* A.strides[0], strides[1]* A.strides[1])+ A.strides
    return as_strided(A, shape= shape, strides= strides)


def block_view_image(A, block=(3, 3), strides=(2, 2)):
    """Provide a 2D tiled block view to a 3d image array. No error checking
    is made.  Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved ;-)
    shape = (A.shape[0] / block[0], A.shape[1] / block[1]) + block
    strides = (
              strides[0] * A.strides[0], strides[1] * A.strides[1]) + A.strides
    return as_strided(A, shape=shape, strides=strides)
