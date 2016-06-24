from numpy.lib.stride_tricks import as_strided
import itertools

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
    strides = (strides[0] * A.strides[0],
               strides[1] * A.strides[1]) + A.strides
    return as_strided(A, shape=shape, strides=strides)

def split_with_overlap(iterable, length, overlap=0, partial=True):
    """
    Return the iterable (array) split to windows of length with a requested
    overlap.  partial controls how the last split is handled.  The options
    are:

    partial=True  - Return any partial array in it's own split
    partial=False - Return any partial array by extending the last full split
    partial=None  - Drop any partial array

    Loosely based off of:
    https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch19s08.html

    Parameters
    ----------
    iterable
        The array or iterable to be split
    length
        length of the splits
    overlap
        overlap to be kept between the splits
    partial
        how to handle partial splits

    Returns
    -------
    list of lists
        A list containing the requested potentially overlapping splits

    """
    it = iter(iterable)

    if partial == True or partial == None:
        results = list(itertools.islice(it, length))
        while len(results) == length:
            yield results
            results = results[length-overlap:]
            results.extend(itertools.islice(it, length-overlap))
        if partial == None:
            pass
        if partial == True:
            if results:
                yield results
    elif partial == False:
        results = list(itertools.islice(it, length*2))
        while True:
            if len(results) < length*2:
                yield results
                break
            else:
                yield results[:length]
                results = results[length-overlap:]
                results.extend(itertools.islice(it, length-overlap))
