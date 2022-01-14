# pinknoise --- a pink-noise wavelet toolbox
----------------------------------------------------------------

`pinknoise` is a library containing a set of tools that implement some of the wavelet-based algorithms tailored to pink-noise (a.k.a., flicker, 1/f-noise) presented in [Carter & Winn (2009)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0909.0747).

**Author**: Nestor Espinoza (nespinoza@stsci.edu)

## Statement of need
Have you ever wondered how to implement likelihood-based inference in the presence of 1/f noise on one-dimensional time-series? This is your library.

## Using the library
Using the library on a given time-series is extremely easy. Let's assume you have a time-series `data`, which hopefully is uniformly sampled. To compute the 
log-likelihood of this dataset given parameters `sigma_w` and `sigma_r`, you would do:

        import pinknoise

        model = pinknoise.compute(len(data))
        model.get_likelihood(data, sigma_w, sigma_r)

Alternatively, you can also define the exponent of a 1/f^gamma pink-noise model:

        model.get_likelihood(data, sigma_w, sigma_r, gamma = 1.1)
        
As per Carter & Winn (2009)'s suggestion, though, don't let the value of gamma go too far off 1.

## Installation
Installation is as simple as:

        python setup.py install

Note this will compile some C code in the background. The basis of the whole wavelet-calculation is done in C.

## Licence and attribution

Read the `LICENCE` file for licencing details on how to use the code. If you make use of this code, please cite 
[Carter & Winn (2009)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0909.0747), and link back to to this repository.
