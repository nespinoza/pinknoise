from .utils import *

class compute(object):
    """
    This class perform the main computations for our time-series.

    :param ndatapoints: (int)
        Number of datapoints in the time-series. This is used to zero-pad the time-series. Note the time-series is expected to have constant sampling. 
    """
    def normal_like(self, x, mu, tau):
        return 0.5*(np.log(tau) - log2pi - tau*( (x-mu)**2))

    def get_fn_likelihood(self, residuals, sigma_w, sigma_r, gamma=1.0):

        like=0.0

        # Generate zero-padded residuals:
        self.zero_padded_residuals[:self.ndatapoints] = residuals

        # Perform wavelet transform of residuals:
        aa, bb = PerformWaveletTransform(self.zero_padded_residuals, self.nzeropadded, self.C, self.nC, self.M)

        # Calculate the g(gamma) factor used in Carter & Winn...
        if(gamma==1.0):
           g_gamma=1.0/(2.0*np.log(2.0))  # (value assuming gamma=1)
        else:
           g_gamma=(2.0)-(2.0)**gamma
        # log-Likelihood of the aproximation coefficients
        sigmasq_S=(sigma_r**2)*g_gamma+(sigma_w)**2
        tau_a =  1.0/sigmasq_S
        like += normal_like( bb[0], 0.0 , tau_a )
        k=long(0)
        SS=range(self.M)
        for ii in SS:
                # log-Likelihood of the detail coefficients with m=i...
                if(ii==0):
                  sigmasq_W=(sigma_r**2)*(2.0**(-gamma*np.double(1.0)))+(sigma_w)**2
                  tau=1.0/sigmasq_W
                  like += normal_like( bb[1], 0.0, tau )
                else:
                  sigmasq_W=(sigma_r**2)*(2.0**(-gamma*np.double(ii+1)))+(sigma_w)**2
                  tau=1.0/sigmasq_W
                  for j in range(2**ii):
                      like += normal_like( aa[k], 0.0 , tau )
                      k=k+1
        return like

    def __init__(self, ndatapoints, wavelet = 'daub4'):

        # Use number of datapoints to create mock dataset:
        self.ndatapoints = ndatapoints
        self.zero_padded_residuals, self.M = zero_pad(self.ndatapoints)
        self.nzeropadded = len(self.zero_padded_residuals)

        # Set the wavelet being used:
        if (wavelet == 'daub4'):
            self.wavelet = 'daub4'
            self.C = np.double(np.arange(4))
            self.C[0] = abs((1.0+sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[1] = abs((3.0+sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[2] = abs((3.0-sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[3] = -abs((1.0-sqrt(3.0))/(4.0*sqrt(2.0)))
            self.nC = 4


