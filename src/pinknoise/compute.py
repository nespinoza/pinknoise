import numpy as np
from .utils import *

log2pi = np.log(2.*np.pi)
class compute(object):
    """
    This class perform the main computations for our time-series.

    :param ndatapoints: (int)
        Number of datapoints in the time-series. This is used to zero-pad the time-series. Note the time-series is expected to have constant sampling. 
    """
    def normal_like(self, x, mu, tau):
        """
        Normal log-likelihood:
        """
        return 0.5*(np.log(tau) - log2pi - tau*( (x-mu)**2))

    def get_likelihood(self, residuals, sigma_w, sigma_r, gamma=1.0):

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

        # Now we calculate the log-Likelihood using equation (32) in Carter & Winn. First, the sigma given by eq. (34):
        sigmasq_S=(sigma_r**2)*g_gamma+(sigma_w)**2
        tau_a =  1.0/sigmasq_S
    
        # The trick is that the FWT code sets n0 in the notation of Carter & Winn to 1 (i.e., the code only works if we have 
        # datasets of size N = n0 * 2**M). So, with this, we only have to compute one term for the second product in eq. (32) 
        # (the one that contains sigma_S^2, sigmasq_S above):
        like += self.normal_like( bb[0], 0.0 , tau_a )

        # Now we have to handle the first term in the product. We compute it sequentially (TODO: this probably can be moved 
        # to a C function to speed up the likelihood evaluation). We iterate over all the possible m's in equation (32); 
        # note ii is *not* m:
        k = 0
        SS=range(self.M)
        for ii in SS:
                # We first do the m = 2. Here n goes from 1 to n0*2 = 2.  
                if(ii==0):
                  sigmasq_W=(sigma_r**2)*(2.0**(-gamma*np.double(1.0)))+(sigma_w)**2
                  tau=1.0/sigmasq_W
                  like += self.normal_like( bb[1], 0.0, tau )
                else:
                  sigmasq_W=(sigma_r**2)*(2.0**(-gamma*np.double(ii+1)))+(sigma_w)**2
                  tau=1.0/sigmasq_W
                  for j in range(2**ii):
                      like += self.normal_like( aa[k], 0.0 , tau )
                      k=k+1
        return like

    def __init__(self, ndatapoints, wavelet = 'daub4'):

        # Use number of datapoints to create mock dataset:
        self.ndatapoints = ndatapoints
        self.zero_padded_residuals, self.M = zero_pad(self.ndatapoints)
        self.nzeropadded = len(self.zero_padded_residuals)

        # Set the wavelet being used:
        if (wavelet == 'daub4'):
            # Coefficients are taken from eq. 13.10.5 in Numerical Recipes:
            self.wavelet = 'daub4'
            self.C = np.double(np.arange(4))
            self.C[0] = abs((1.0+sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[1] = abs((3.0+sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[2] = abs((3.0-sqrt(3.0))/(4.0*sqrt(2.0)))
            self.C[3] = abs((1.0-sqrt(3.0))/(4.0*sqrt(2.0)))
            self.nC = 4


