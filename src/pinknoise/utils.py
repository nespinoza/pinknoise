import numpy as np
from math import sqrt,log
import FWT

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
