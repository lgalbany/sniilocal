
# coding: utf-8

# In[1]:

#!/usr/bin/env python


"""This program pre-processes host galaxy spectra for use in STARLIGHT.
The data format (MANGA OR PMAS) must be specified for the proper corrections
to be made to the data. The data is read in from the .fits file, and
depending on the data format, corrections are made to account for 
Milky Way extinction, flux reddening, and redshift. This code produces
the correctly formatted input files for STARLIGHT so after a .fits file
is processed through this program, STARLIGHT may be run.
"""

from collections import namedtuple
from astropy.io import fits
from astropy.table import Table
from extinction import fitzpatrick99
import os
import astropy
import extinction
import numpy as np
import argparse


__author__ = 'Meghan Cilento'


# In[2]:

def generate_run_parameters(survey_name, path):
    """Indicate the set of parameters specific for the data format of the .fits
    file being run through this program.
    
    This program runs with two different file formats, MANGA and PMAS. This
    must be specified by the user. Once the file format is specified, the
    program will run with the parameters set in this function.
    
    Args:
        survey_name : the name of the survey conducted (MANGA or PMAS)
        path        : the name of the fits file
                      (this should be in the same directory as this program)
        
    Returns:
        parameters(): the set of parameters dependent on the survey_name
    """
    
    with fits.open(path) as hdulist:
        parameters = namedtuple('Parameters', ['wavelength_axis', 'step', 'initial_wavelength', 'flux_index', 'eflux_index', 'r_v', 'factor'])
        
        if survey_name == "MANGA":
            wavelength_axis = hdulist[1].header['NAXIS3']
            step = hdulist[1].header['CD3_3']
            initial_wavelength = hdulist[1].header['CRVAL3']
            
            return parameters(wavelength_axis=wavelength_axis, step=step, initial_wavelength=initial_wavelength, flux_index=1, eflux_index=2, r_v=3.1, factor=1.0e-17)

        elif survey_name == "PMAS":
            wavelength_axis = hdulist[0].header['NAXIS3']
            step = hdulist[0].header['CDELT3']
            initial_wavelength = hdulist[0].header['CRVAL3']
            
            return parameters(wavelength_axis=wavelength_axis, step=step, initial_wavelength=initial_wavelength, flux_index=0, eflux_index=1, r_v=3.1, factor=1.0e-16)


# In[3]:

def generate_resampled_data(HDULIST, z, a_v, params, survey_name):
    """Generate the wavelengths, flux, and flux error used in STARLIGHT input
    files.
    
    This function reads in wavelength, flux, and flux error data from the fits
    file. The wavelengths are shifted into rest frame and resampled uniformly
    with a step of 1 angstrom. The flux and flux error data is shifted into
    rest frame, dereddened, and resampled uniformly. The resampled arrays will
    appear in the STARLIGHT input files in the form of columns.
    
    Args:
        HDULIST : the opened fits file from which data is read in
        z       : the redshift parameter
        a_v     : extinction parameter specific for each data set
        r_v     : Milky Way extinction parameter - fixed value of 3.1
        
    Returns:
        resampled_wavelength  : array of wavelengths resampled
        resampled_flux        : array of flux values resampled
        resampled_eflux       : array of flux error values resampled
    """
    
    last_wavelength = params.initial_wavelength + (params.step * params.wavelength_axis)
    num_of_elements = ((last_wavelength - params.initial_wavelength) / params.step)
    wavelength = np.linspace(params.initial_wavelength, last_wavelength, num_of_elements, endpoint=True)
    
    restframed_wavelength = wavelength / (1 + z)
    number_of_wavelengths = np.floor(max(restframed_wavelength)) - np.ceil(min(restframed_wavelength))
    resampled_wavelength = np.arange(number_of_wavelengths) + np.ceil(restframed_wavelength[0])

    
    flux = HDULIST[params.flux_index].data * params.factor
    eflux = np.sqrt(1./ HDULIST[params.eflux_index].data) * params.factor
    
    if survey_name == "MANGA":
        a_lambda = extinction.fitzpatrick99(wavelength, a_v, params.r_v)
        correction_factor = 10**(0.4*a_lambda)
        
    elif survey_name == "PMAS":
        correction_factor = 1
        
    resampled_flux = np.empty((len(resampled_wavelength), flux.shape[1], flux.shape[2]))
    resampled_eflux = np.empty((len(resampled_wavelength), flux.shape[1], flux.shape[2]))
    
    for i in range(flux.shape[1]):
        for j in range(flux.shape[2]):
            unreddened_flux = flux[:, i, j] * (1 + z) * correction_factor
            unreddened_eflux = eflux[:, i, j] * (1 + z) * correction_factor
            
            resampled_flux[:, i, j] = np.interp(resampled_wavelength, 
                                                restframed_wavelength, 
                                                unreddened_flux)
            resampled_eflux[:, i, j] = np.interp(resampled_wavelength, 
                                                restframed_wavelength,
                                                unreddened_eflux)
            
    return resampled_wavelength, resampled_flux, resampled_eflux


# In[4]:

def generate_flag(z, resampled_wavelength):
    """Generate the flag array used in STARLIGHT input files.
    
    This function creates a flag array used to mask irrelevant data points
    in the flux and flux error arrays. The wavelengths are fixed for any 
    data set run through this program. Do not change them. Value of 3 masks
    the data point, value of 0 uses the data point.
    
    Arguments:
        z              : the redshift parameter
        resampled_wavelength : array of wavelengths resampled
        
    Returns:
        flag : fixed array which masks sections of the input spectra
    """
    
    bad = np.array([5453, 5468,
      5565, 5589,
      5883, 5903,
      6294, 6306,
      4350, 4366,
      7162, 7333, 
      9900, 11000]) 

    bad_restframed = bad / (1 + z)
    a = bad_restframed.reshape(7, 2)
    
    flag = np.zeros(len(resampled_wavelength), dtype = np.int8)

    for i in range(0, 7):
        flag[(a[i][0] < resampled_wavelength) & (a[i][1] > resampled_wavelength)] = 3

    return flag


# In[5]:

def generate_condition_lists(resampled_wavelength, resampled_flux):
    """Generate the signal to noise values and corresponding indices for each
    pixel.
    
    This function creates two arrays, one filled with the indices
    corresponding to signal to noise values > 5, the other filled 
    with the pixel indices in which 95% of the resampled flux array has
    values greater than 0. If less than 95% is greater than 0, then the pixel
    is not used. This function is necessary to filter through pixels in 
    order to use the correct ones. Only pixels with signal to noise > 5
    and 95% flux > 0 are relevant. 
    
    Arguments:
        resampled_wavelength : array of wavelengths resampled
        resampled_flux       : array of flux values resampled
        
    Returns:
        S_N_indices : list of all indices where signal to noise > 5
        good_pixels : list of pixel indices where 95% of resampled flux > 0
    """
    
    S_N_values = []
    S_N_indices = []
    good_pixels = []
    
    for i in range(resampled_flux.shape[1]):
        for j in range(resampled_flux.shape[2]):
            
            resampled_flux_i_j = resampled_flux[1091:1151, i, j]

            coefficent = np.polyfit(resampled_wavelength[1091:1151], resampled_flux_i_j, 3)

            x = resampled_wavelength[1091:1151]

            equation = (coefficent[0]*(x**3)) + (coefficent[1]*(x**2)) + (coefficent[2]*x) + (coefficent[3])
            
            S_N = 1. / np.std((resampled_flux_i_j / equation), ddof = 1)
            
            S_N_values.append(S_N)
            
            if S_N > 5:
                S_N_indices.append('%s_%s'%(i, j))
                
            if np.count_nonzero(resampled_flux[:, i, j]) > (0.95*len(resampled_wavelength)):
                good_pixels.append('%s_%s'%(i, j))

    return S_N_indices, good_pixels


# In[6]:

def generate_header_info(resampled_wavelength, S_N_indices, good_pixels, filename):
    """Generate the header information used in the STARLIGHT input file.
    
    This function creates the top 15 lines of a STARLIGHT input file. This file
    will be called in generate_STARLIGHT_input(), which writes the input file
    used in STARLIGHT. Every input file must have these first 15 lines.
    
    Arguments:
        resampled_wavelength : array of wavelengths resampled
        S_N_indices          : list of all indices where signal to noise > 5
        good_pixels          : list of pixel indices where 95% of resampled
                               flux > 0

    Returns:
        header               : first 15 lines included at the top of each grid
                               file
    """
    
    # Arbitrarily choose +/- 50 to include the fringes at the end of the
    # spectrum
    fits_num = len(set(S_N_indices) & set(good_pixels))
    min_wave = min(resampled_wavelength) + 50
    #max_wave = max(resampled_wavelength) - 50 (instead default to 7500)
    
    header = "{}                                    [Number of fits to run]\n"\
             "./modelsnew/                           [base_dir]\n"\
             "./{} files                                     [obs_dir]\n"\
             "./                                     [mask_dir]\n"\
             "./{} files                                     [out_dir]\n"\
             "-2007200                               [your phone number]\n"\
             "4580.0                                 [llow_SN]   lower-lambda of S/N window\n"\
             "4640.0                                 [lupp_SN]   upper-lambda of S/N window\n"\
             "{}                                 [Olsyn_ini] lower-lambda for fit\n"\
             "7500                                 [Olsyn_fin] upper-lambda for fit\n"\
             "1.0                                    [Odlsyn]    delta-lambda for fit\n"\
             "1.0                                    [fscale_chi2] fudge-factor for chi2\n"\
             "FIT                                    [FIT/FXK] Fit or Fix kinematics\n"\
             "1                                      [IsErrSpecAvailable] 1/0 = Yes/No\n"\
             "1                                      [IsFlagSpecAvailable] 1/0 = Yes/No\n"\
        .format(fits_num, filename, filename, min_wave)
    
    return header


# In[7]:

def generate_STARLIGHT_input(resampled_wavelength, resampled_flux, resampled_eflux, flag, header, filename):
    """Generate the STARLIGHT input file.
    
    This function writes the STARLIGHT input file. It contains the 15 lines
    created with the generate_header_info() function and columns of text files,
    each of which correspond to pixels that have 95% of flux greater than 0 and
    signal to noise greater than 5. Each file contains wavelength, flux, error,
    and flag columns corresponding to the spectrum of that particular pixel.
    
    Arguments:
        resampled_wavelength : array of wavelengths resampled
        resampled_flux       : array of flux values resampled
        resampled_eflux      : array of flux error values resampled
        flag                 : fixed array which masks sections of the input
                               spectra
        header               : first 15 lines included at the top of each grid
                               file
        filename             : name of the fits file (only beginning)
    """
    
    if not os.path.exists('./%s files'%(filename)):
        os.mkdir('./%s files'%(filename))
        
    table = Table(names=['file_paths', 'StCv04.C11.config', 'Base.CB07.vall.4m', 
                         'Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0', 'file_paths_again'],
                  dtype=[object, object, object, object, object, object, object, object]) 
    
    for i in range(resampled_flux.shape[1]):
        for j in range(resampled_flux.shape[2]):

            resampled_flux_i_j = resampled_flux[:, i, j]

            if np.count_nonzero(resampled_flux_i_j) > (0.95*len(resampled_wavelength)):

                resampled_flux_i_j = resampled_flux[1091:1151, i, j]

                coefficent = np.polyfit(resampled_wavelength[1091:1151], resampled_flux_i_j, 3)

                x = resampled_wavelength[1091:1151]

                equation = (coefficent[0]*(x**3)) + (coefficent[1]*(x**2)) + (coefficent[2]*x) + (coefficent[3])

                S_N = 1. / np.std((resampled_flux_i_j / equation), ddof = 1)

                if S_N > 5:
                    resampled_flux_i_j = resampled_flux[:, i, j]
                    resampled_eflux_i_j = resampled_eflux[:, i, j]
                    t = Table([resampled_wavelength, resampled_flux_i_j, resampled_eflux_i_j, flag], 
                              names = ['wavelength', 'flux', 'error', 'flag'])
                    t.write('./%s files/%s_%s_%s.txt'%(filename, filename, i, j),
                            format = 'ascii', overwrite=True)

                    table.add_row(['%s_%s_%s.txt'%(filename, i, j), 'StCv04.C11.config', 
                                   'Base.CB07.vall.4m', 'Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0', 
                                   '%s_%s_%s.C11.gm.CCM.BN'%(filename, i, j)])

    with open('./%s_grid.in'%(filename), 'w') as outfile:
        outfile.write(header)

        for row in table:
            outfile.write(' '.join(row) + '\n')


# In[8]:

def pre_process_image_data(z, a_v, survey_name, filename, path):
    """Pre-process a fits file for use in STARLIGHT.
    
    This function completely processes a fits file whose spectra will be used
    in STARLIGHT. All other functions are called and the file output of this
    function is the input grid file, which is formatted by the standards set in
    the STARLIGHT manual. Once this program runs, the grid file ready to
    be run through STARLIGHT.

    Arguments:
        z           :   the redshift parameter
        a_v         :   extinction parameter specific for each data set
        survey_name :   the name of the survey conducted (MANGA or PMAS)
        filename    :   name of the fits file (only beginning)
        path        :   the name of the fits file
                        (this should be in the same directory as this program)
    """
    
    HDULIST = fits.open(path)
    
    params = generate_run_parameters(survey_name, path)

    resampled_wavelength, resampled_flux, resampled_eflux = generate_resampled_data(HDULIST, z, a_v, params, survey_name)
    
    flag = generate_flag(z, resampled_wavelength)
    
    S_N_indices, good_pixels = generate_condition_lists(resampled_wavelength, resampled_flux)
    
    header = generate_header_info(resampled_wavelength, S_N_indices, good_pixels, filename)
    
    generate_STARLIGHT_input(resampled_wavelength, resampled_flux, resampled_eflux, flag, header, filename)



# In[ ]:

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-z', '--redshift', type=float, default=None)
PARSER.add_argument('-a', '--a_v', type=float, default=None)
PARSER.add_argument('-n', '--survey_name', type=str, default=None)
PARSER.add_argument('-f', '--filename', type=str, default=None)
PARSER.add_argument('-p', '--path', type=str, default=None)

if __name__ == '__main__':
    args = PARSER.parse_args()
    pre_process_image_data(args.redshift, args.a_v, args.survey_name, args.filename, args.path)


