def read_IFU_cube_noERR(cube, survey, path_cubes):
    
    ''' 
        Read an IFU cube and process its data
        Surveys: PMAS/MUSE
        Stellar populations (starsub) are either subtracted or not
        Remember: $ heainit  -->  $fv cube.fits & exit  -->  Explore cube from the terminal
        Example: NGC2906
        
        ------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------
        
        PMAS:: shape: (3773, 71, 78) --> Transpose to original: (78, 71, 3773)
        PMAS[0] = DATA    -->  Flux
        PMAS[1] = ERROR   -->  Error
        PMAS[2] = ... [many 0's and some 1's]
        
        PMAS[0].header['CRVAL3'] : starting wavelength          = 3675 A
        PMAS[0].header['CDELT3']/['CD3_3'] : wavelength step    = 1 A
        PMAS[0].header['NAXIS3'] : number of wavelength steps   = 3773
        
        ------------------------------------------------------------------------------------
        
        starsub_PMAS:: shape: (3773, 71, 78) --> Transpose to original: (78, 71, 3773)
        PMAS[0] = DATA    -->  Flux
        PMAS[1] = ERROR   -->  Error
        
        starsub_PMAS[0].header['CRVAL3'] : starting wavelength          = 3675 A
        starsub_PMAS[0].header['CDELT3']/['CD3_3'] : wavelength step    = 1 A
        starsub_PMAS[0].header['NAXIS3'] : number of wavelength steps   = 3773
        
        ------------------------------------------------------------------------------------
        
        MUSE:: shape: (4567, 322, 320) --> Transpose to original: (320, 322, 4567)
        MUSE[0] = DATA    -->  Flux
        MUSE[1] = ERROR   -->  Error (variance)
        
        MUSE[1].header['CRVAL3'] : starting wavelength          = 4717 A
        MUSE[1].header['CDELT3']/['CD3_3'] : wavelength step    = 1 A
        MUSE[1].header['NAXIS3'] : number of wavelength steps   = 3681
        
        ------------------------------------------------------------------------------------
        
        starsub_MUSE:: shape: (4567, 322, 320) --> Transpose to original: (320, 322, 4567)
        starsub_MUSE[0] = DATA    -->  Flux
        starsub_MUSE[1] = ERROR   -->  Error (variance)
        
        starsub_MUSE[0].header['CRVAL3'] : starting wavelength          = 4717 A
        starsub_MUSE[0].header['CDELT3']/['CD3_3'] : wavelength step    = 1 A
        starsub_MUSE[0].header['NAXIS3'] : number of wavelength steps   = 4567
    
    '''
    
    import numpy as np
    import os
    import astropy.io.fits as fits
    
    import warnings
    warnings.filterwarnings("ignore") # Ignore astropy.io.fits warning with starsub PMAS:
    # WARNING: The following header keyword is invalid or follows an unrecognized non-standard convention:
    # ADDED MANUALLY MED_VEL = ...
    
    
    
    #####################################################################################################
    
    
    
    if survey not in ('PMAS', 'MUSE'): # Add other surveys if desired
        
        raise NameError('Please, analyze a valid survey: PMAS/MUSE')
        
    else: # Main code begins


    
        readcube = fits.open(path_cubes + cube + '.fits', format = 'fits')


        # Copy header for MPFIT cube
        hdr = readcube[0].header


        # Avoid some problems with headers
        if 'COMMENT' in hdr:

            hdr.remove('COMMENT', remove_all = True)

    

    
        ########## Decompose the cube to retrieve wavelength, flux and error for each spaxel ##########

        wv_start = hdr['CRVAL3']
        wv_nelem = hdr['NAXIS3']

    
        if ( ('CDELT3' in hdr) | (hdr['CD3_3'] == 0.) ):

            wv_step =  hdr['CDELT3']


        elif ( ('CD3_3' in hdr) | (hdr['CDELT3'] == 0.) ):

            wv_step =  hdr['CD3_3']

        else:
        
            raise NameError('Please, analyze a valid survey: PMAS/MUSE')




        wavelength = np.asarray([wv_start + wv_step*x for x in range(wv_nelem)])




        Flux = readcube[0].data.T
        #Error = readcube[1].data.T
    
    
    
    
        # Corrections based on a) surveys and b) starsub
        # Flux and error will have units of [ergs/cm2/angstrom/sec]. Typical values ~10^{-16} ergs/cm2/angstrom/sec
    
        if ( (survey == 'PMAS') & ('resam' in cube) ): # PMAS, resam
        
            Flux = 1.E-16 * Flux                        #[1E-16 ergs/cm2/angstrom/sec]
            #Error = 1.E-16 * Error                      #[1E-16 ergs/cm2/angstrom/sec]
        
        elif ( (survey == 'PMAS') & ('starsub' in cube) ): # PMAS, starsub
        
            pass #[ergs/cm2/angstrom/sec]
        
        elif ( (survey == 'MUSE') & ('resam' in cube) ): # MUSE, resam
        
            Flux = 1.E-20 * Flux                        #[1E-20 ergs/cm2/angstrom/sec]
            #Error = 1.E-20 * np.sqrt(Error)             #[1E-20 ergs/cm2/angstrom/sec]**2 (VARIANCE) -> SIGMA = NP.SQRT(VARIANCE)
        
        elif ( (survey == 'MUSE') & ('starsub' in cube) ): # MUSE, starsub
        
            #Error = 1.E-20 * np.sqrt(1.E20 * Error)     # Another conversion needed... 
	    pass

        else: # star_norm, equivalent widths

            pass 




        # Correct nans
        wh_nan = np.where(np.isfinite(Flux) == False)
        Flux[wh_nan] = 1.E-22
        #wh_nan = np.where(np.isfinite(Error) == False)
        #Error[wh_nan] = 1.E-22


        #return [wavelength, Flux, Error, hdr]
        return [wavelength, Flux, hdr]
