def cut_and_normalize_spectrum(w, flux, limits):
    '''Normalize spectrum between 0 and 1 in a desired wavelength range
    Originally thought for CALIFA/MUSE IFU spectra

    Parameters
    ----------
    w : 1 dim np.ndarray
    array of wavelengths
    flux : np.ndarray of shape (len(w))
    array of flux values
    limits : list
    wavelengths for final spectrum [low1,up1]
    Example:
    wHalpha, fHalpha = cut_and_normalize_spectrum(wavelength, flux,
    [[6500*u.AA, 6620*u.AA])
    '''
    import numpy as np
    # useful when dealing with columns
    w = np.asarray(w)
    flux = np.asarray(flux)
    #index of the region we want to return
    indrange = (w > limits[0]) & (w < limits[1])
    # make a flux array of shape
    # (number of spectra, number of points in indrange)
    f = np.zeros((indrange.sum()))
    # subtract the local continuum
    f = (flux[indrange] - np.min(flux[indrange]))/np.max((flux[indrange] - np.min(flux[indrange])))
    return w[indrange], f




def cut_and_normalize_spectrum_cc(w, flux, limits):
    '''Normalize spectrum between 0 and 1 in a desired wavelength range
    Originally thought for CALIFA/MUSE IFU spectra
    LLUIS' CODE FOR THE CROSS-CORRELATION COEFFICIENTS

    Parameters
    ----------
    w : 1 dim np.ndarray
    array of wavelengths
    flux : np.ndarray of shape (len(w))
    array of flux values
    limits : list
    wavelengths for final spectrum [low1,up1]
    Example:
    wHalpha, fHalpha = cut_and_normalize_spectrum(wavelength, flux,
    [[6500*u.AA, 6620*u.AA])
    '''
    import numpy as np
    # useful when dealing with columns
    w = np.asarray(w)
    flux = np.asarray(flux)
    #index of the region we want to return
    indrange = (w > limits[0]) & (w < limits[1])
    # make a flux array of shape
    # (number of spectra, number of points in indrange)
    f = np.zeros((indrange.sum()))
    # normalize
    #thres = np.mean(flux[indrange][0:21])
    f = (flux[indrange] - np.min(flux[indrange]))/np.sum((flux[indrange] - np.min(flux[indrange])))
    #f = (flux[indrange] - thres)/np.sum((flux[indrange] - thres))
    return w[indrange], f