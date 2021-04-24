import numpy as np
from math import sqrt,log

from stochastic.processes.noise import ColoredNoise

import FWT

def generate_detectorTS(columns, rows, pixel_time = 10, jump_time = 120, return_image = False, beta = 1., sigma = 1., add_poisson = False, poisson_sigma = 1.):
    """
    This function simulates a JWST detector image and corresponding time-series of the pixel-reads, assuming the noise follows a $1/f^\beta$ power-law in its 
    power spectrum. This assumes the 1/f pattern (and hence the detector reads) go along the columns of the detector.

    Parameters
    ----------
    columns : int
        Number of columns of the detector.

    rows : int
        Number of rows of the detector.

    pixel_time : float
        Time it takes to read a pixel along each column in microseconds. Default is 10 microseconds (i.e., like JWST NIR detectors).

    jump_time : float
        Time it takes to jump from one column to the next once all its pixels have been read, in microseconds. Default is 120 microseconds (i.e., like JWST NIR detectors).

    return_image : boolean
        If True, returns an image with the simulated values. Default is False.

    beta : float
        Power-law index of the PSD of the noise. Default is 1.

    sigma : float
        Variance of the power-law process in the time-domain. Default is 1.

    add_poisson : boolean
        If True, poisson noise with mean and variance equal to the square of `poisson_sigma` will be (additively) added to the power-law process. Default is False.

    poisson_sigma : boolean
        Square-root of the variance (and hence of the mean) of the added Poisson noise if `add_poisson` is True. Default is 1.

    Returns
    -------
    times : `numpy.array`
        The time-stamp of the flux values (i.e., at what time since read-out started were they read).

    time_series : `numpy.array`
        The actual flux values on each time-stamp (i.e., the pixel counts as they were read in time).

    image : `numpy.array` 
        The image corresponding to the `times` and `time_series`, if `return_image` is set to True.

    """
    # This is the number of "fake pixels" not read during the waiting time between jumps:
    nfake = int(jump_time/pixel_time)

    # First, generate a time series assuming uniform sampling (we will chop it later to accomodate the jump_time):
    CN = ColoredNoise(beta = beta, t = (rows * columns * pixel_time) + columns * jump_time)

    # Get the samples and time-indexes:
    nsamples = rows * columns + (nfake * columns)
    y = CN.sample(nsamples)
    t = CN.times(nsamples)

    # Now remove samples not actually read by the detector due to the wait times:
    all_indexes = np.arange(nsamples)
    idx = []
    for i in range(columns):

        min_idx = i * (rows + nfake)
        idx = idx + list(all_indexes[min_idx:min_idx + rows])

    times = t[idx]
    time_series = y[idx]

    # Set process standard-deviation to input sigma:
    time_series = sigma * (time_series / np.sqrt(np.var(time_series)) )

    # Add poisson noise if user wants it:
    if add_poisson:
        time_series = time_series + np.random.poisson(poisson_sigma**2, len(time_series))

    if return_image:
        # Create image:
        image = time_series.reshape((columns, rows)).transpose()
        # Return all:
        return times, time_series, image
    else:
        return times, time_series

def zero_pad(ndata, data = None):
    """
    This function creates a mock array of length ndata + nextra optimized for fast wavelength transform 
    calculations. The nextra elements at the end are zeroes. Function returns this mock array, 
    and M --- the scale of the wavelet (in base 2).
    """

    fdatalen=0
    M=0

    # First we search for the optimal data length...
    for i in range(ndata):

        min_2 = 2**i
        max_2 = 2**(i+1)
        if( min_2 < ndata and max_2 > ndata ):
            fdatalen = max_2
            M = i + 1

        elif( min_2 == ndata ):
            fdatalen = ndata
            M = i

        if( fdatalen!=0 ):
            break

    # Now we form our zero-padded vector (we padd the zeroes to the end of the vector)...
    mock_data_vector = np.zeros(fdatalen)
    if data is None:
        mock_data_vector[:ndata] = 1.
    else:
        mock_data_vector[:ndata] = np.copy(data)

    # Return mock vector and M
    return mock_data_vector, M

def PerformWaveletTransform(data_vector, ndatavector, C, nC, M):
    """
    This function performs the wavelet transform of a vector of length 2**i, where i is any integer. 
    If your dataset does not have this length, pad zeroes at the end with the zero_pad function. 
    """

    Result1, Result2 = FWT.getWC(data_vector, C, ndatavector, nC, M)
    FinalMatrix1 = np.asarray(Result1)
    FinalMatrix2 = np.asarray(Result2)

    return FinalMatrix1, FinalMatrix2
 
def PerformInverseWaveletTransform(data_vector, ndatavector, C, nC, M):
    """
    This functions performs the inverse wavelength transform of a data vector
    """

    Result = FWT.getSignal(data_vector, C, ndatavector, nC, M)
    Signal = np.asarray(Result)

    return Signal
