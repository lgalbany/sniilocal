def basic_read_IFU_cube(cube):
    
    ''' Read an IFU cube without processing its data'''

    import os
    import astropy.io.fits as fits
    
    import warnings
    warnings.filterwarnings("ignore") # Ignore astropy.io.fits warning with starsub CALIFA:
    # WARNING: The following header keyword is invalid or follows an unrecognized non-standard convention:
    # ADDED MANUALLY MED_VEL = ...

    path_cubes = os.getcwd() + '/Cubes/'
    
    readcube = fits.open(path_cubes + cube + '.fits', format = 'fits')
    
    return readcube









def read_IFU_cube_CALIFAMUSE(cube, survey):#, starsub = False):
    
    ''' 
    Read an IFU cube and process its data
    Surveys: CALIFA/MUSE
    Stellar populations (starsub) are either subtracted or not
    Remember: $ heainit  -->  $fv cube.fits & exit  -->  Explore cube from the terminal
    
    ------------------------------------------------------------------------------------
    ------------------------------------------------------------------------------------
    
    CALIFA:: shape: (1901, 71, 78) --> Transpose to original: (78, 71, 1901)
    CALIFA[0] = Primary      -->  Header, flux
    CALIFA[1] = ERROR        -->  Error
    CALIFA[2] = ERRWEIGHT    --> Weighted error
    CALIFA[3] = BADPIX       --> Bad pixels
    CALIFA[4] = FIBCOVER     --> Fiber coverage
    
    CALIFA[0].header['CRVAL3'] : starting wavelength          = 3701 A
    CALIFA[0].header['CDELT3'] : wavelength step              = 2 A
    CALIFA[0].header['NAXIS3'] : number of wavelength steps   = 1901
    
    ------------------------------------------------------------------------------------
    
    starsub_CALIFA:: shape: (3773, 71, 78) --> Transpose to original: (78, 71, 3773)
    CALIFA[0] = Primary  -->  Flux
    CALIFA[1] = NoName   -->  Error (variance)
    
    starsub_CALIFA[0].header['CRVAL3'] : starting wavelength          = 3675 A
    starsub_CALIFA[0].header['CDELT3'] : wavelength step              = 1 A
    starsub_CALIFA[0].header['NAXIS3'] : number of wavelength steps   = 3773
    
    ------------------------------------------------------------------------------------
    
    MUSE:: shape: (3681, 322, 320) --> Transpose to original: (320, 322, 3681)
    MUSE[0] = Primary  -->  Header
    MUSE[1] = DATA     -->  Flux
    MUSE[2] = STAT     -->  Error (variance)
    
    MUSE[1].header['CRVAL3'] : starting wavelength          = 4750 A
    MUSE[1].header['CD3_3'] : wavelength step               = 1.25 A
    MUSE[1].header['NAXIS3'] : number of wavelength steps   = 3681
    
    ------------------------------------------------------------------------------------
    
    starsub_MUSE:: shape: (4567, 322, 320) --> Transpose to original: (320, 322, 4567)
    starsub_MUSE[0] = Primary  -->  Flux
    starsub_MUSE[1] = NoName   -->  Error (variance)
    
    starsub_MUSE[0].header['CRVAL3'] : starting wavelength          = 4717 A
    starsub_MUSE[0].header['CD3_3'] : wavelength step               = 1.25 A
    starsub_MUSE[0].header['NAXIS3'] : number of wavelength steps   = 4567
    
    '''
    
    import numpy as np
    import os
    import astropy.io.fits as fits
    
    import warnings
    warnings.filterwarnings("ignore") # Ignore astropy.io.fits warning with starsub CALIFA:
    # WARNING: The following header keyword is invalid or follows an unrecognized non-standard convention:
    # ADDED MANUALLY MED_VEL = ...

    path_cubes = os.getcwd() + '/Cubes/'
    
    readcube = fits.open(path_cubes + cube + '.fits', format = 'fits')
    
    
    # Decompose the cube, retrieve wavelength, etc.
    
    if survey == 'CALIFA': #'CALIFA' in cube:
        
        wv_start = readcube[0].header['CRVAL3']
        wv_step =  readcube[0].header['CDELT3']
        wv_nelem = readcube[0].header['NAXIS3']
        wavelength = np.asarray([wv_start + wv_step*x for x in range(wv_nelem)])
        
        if 'starsub' not in cube: #not starsub:
        
            Flux = readcube[0].data.T
            Error = readcube[1].data.T
            Errweight = readcube[2].data.T
            Badpix = readcube[3].data.T
            Fibcover = readcube[4].data.T

            # Correct nans
            wh_nan = np.where(np.isfinite(Flux)==False)
            Flux[wh_nan] = 0.
            wh_nan = np.where(np.isfinite(Error)==False)
            Error[wh_nan] = 0.
            
            #print 'Reading CALIFA cube (unsubtracted stellar populations)...'
            #print '78 columns X 71 rows X 1901 wavelength steps'
            return [wavelength, Flux, Error, Errweight, Badpix, Fibcover]
        



        elif 'starsub' in cube: #starsub:
            
            Flux = readcube[0].data.T
            Error = readcube[1].data.T

            # Correct nans
            wh_nan = np.where(np.isfinite(Flux)==False)
            Flux[wh_nan] = 0.
            wh_nan = np.where(np.isfinite(Error)==False)
            Error[wh_nan] = 0.
        
            #print 'Reading CALIFA cube (subtracted stellar populations)...'
            #print '78 columns X 71 rows X 3773 wavelength steps'
            return [wavelength, Flux, Error]
        
        
        
        
    elif survey == 'MUSE': #'MUSE' in cube:
        
        if 'starsub' not in cube: #not starsub:
            
            wv_start = readcube[1].header['CRVAL3']
            wv_step =  readcube[1].header['CD3_3']
            wv_nelem = readcube[1].header['NAXIS3']
            wavelength = np.asarray([wv_start + wv_step*x for x in range(wv_nelem)])
        
            Flux = readcube[1].data.T
            Error = readcube[2].data.T

            # Correct nans
            wh_nan = np.where(np.isfinite(Flux)==False)
            Flux[wh_nan] = 0.
            wh_nan = np.where(np.isfinite(Error)==False)
            Error[wh_nan] = 0.
        
            #print 'Reading MUSE cube (unsubtracted stellar populations)...'
            #print '320 columns X 322 rows X 3861 wavelength steps'
            return [wavelength, Flux, Error]



        
        elif 'starsub' in cube: #starsub:
            
            wv_start = readcube[0].header['CRVAL3']
            wv_step =  readcube[0].header['CD3_3']
            wv_nelem = readcube[0].header['NAXIS3']
            wavelength = np.asarray([wv_start + wv_step*x for x in range(wv_nelem)])
            
            Flux = readcube[0].data.T
            Error = readcube[1].data.T

            # Correct nans
            wh_nan = np.where(np.isfinite(Flux)==False)
            Flux[wh_nan] = 0.
            wh_nan = np.where(np.isfinite(Error)==False)
            Error[wh_nan] = 0.
        
            #print 'Reading MUSE cube (unsubtracted stellar populations)...'
            #print '320 columns X 322 rows X 4567 wavelength steps'
            return [wavelength, Flux, Error]
        
        
        
        
    else:
        
        raise NameError('Please, introduce a valid survey: CALIFA or MUSE')