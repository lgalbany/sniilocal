
# coding: utf-8

# In[1]:

#!/usr/bin/env python


"""This code processes the output files created by STARLIGHT. There is a
single output file corresponding to each pixel in the host galaxy image. These
pixels are the same as the pixels used in the input files for STARLIGHT. Data
is read in from these files, and from the data, plots of stellar velocity,
stellar age distribution, and spectral fits may be created. These are important
for supernova and host galaxy correlation studies.
"""

import numpy as np
import glob 
from collections import namedtuple
from astropy.io import fits
from  matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

__author__ = 'Meghan Cilento'


# In[2]:

def gen_run_parameters(survey_name, PATH):
    """Indicate the set of parameters specific for the format of the fits file being run
    through this program.
    
    This program runs with two different file formats, MANGA and PMAS. This must be specified
    by the user. Once the file format is specified, the program will run with the parameters
    specifically set in this function. 
    
    Args:
        survey_name : the name of the survey conducted (MANGA or PMAS)
        PATH        : the name of the fits file (this should be in the same directory as this program)
        
    Returns:
        parameters(): the different sets of possible parameters dependent on
                      the survey_name
    """
    
    with fits.open(PATH) as hdulist:
        parameters = namedtuple('Parameters', ['flux_index'])
        
        if survey_name == "MANGA":
            return parameters(flux_index=1)

        elif survey_name == "PMAS":
            return parameters(flux_index=0)


# In[3]:

def generate_data_shape(params, PATH):
    """Generate the dimensions of the data cube.
    
    This function automatically reads the data cube used in the pre-processing program to find 
    the dimensions of the data. This is necessary so that when reading the output files the 
    program knows the indicies through which to iterate and gather pertinent data for analysis.  
    """
    
    HDULIST = fits.open(PATH)
    
    flux = HDULIST[params.flux_index].data
    dim_one = flux.shape[1]
    dim_two = flux.shape[2]
    
    return dim_one, dim_two


# In[4]:

def generate_av_ov_vd_values(dim_one, dim_two):
    """Generate av, ov, and vd values for all pertinent pixels in a host galaxy image.  
    
    This function reads all of the output files and finds the specific values of av, ov, and vd
    within each file. These values are stored in arrays with the same dimensions of the data cube.
    By storing the values with the same dimension, they can be easily plotted with pixel axis and 
    color bars to indicate trends and host galaxy properties.
    """
    
    av_values = np.empty((dim_one, dim_two))
    vo_values = np.empty((dim_one, dim_two))
    vd_values = np.empty((dim_one, dim_two))
    
    av_values[:]=np.nan
    vo_values[:]=np.nan
    vd_values[:]=np.nan
    
    for i in range(dim_one):
        for j in range(dim_two):

            my_files = glob.glob('/manga-8158-1901 files/manga-8158-1901_%s_%s.C11.gm.CCM.BN'%(i, j))

            if my_files:
                for name in my_files:
                    f = open(name, 'r')
                    
                    first_line = f.readline()
                    
                    if len(first_line) == 68:

                        for line in f:

                            if len(line.split()) > 2:
                                if line.split()[2] == '(mag)]' and line.split()[1] == '[AV_min': 
                                    av_values[i, j] = float(line.split()[0])

                            if len(line.split()) > 2:
                                if line.split()[2] == '(km/s)]' and line.split()[1] == '[v0_min': 
                                    vo_values[i, j] = float(line.split()[0])

                            if len(line.split()) > 2:
                                if line.split()[2] == '(km/s)]' and line.split()[1] == '[vd_min': 
                                    vd_values[i, j] = float(line.split()[0])
                                
    return av_values, vo_values, vd_values


# In[26]:

def generate_normalized_light_ages(dim_one, dim_two):

    normalization_values = np.empty((dim_one, dim_two))
    normalization_values[:] = np.nan
    
    for i in range(dim_one):
        for j in range(dim_two):
    
            my_files = glob.glob('./manga-8158-1901 files/manga-8158-1901_%s_%s.C11.gm.CCM.BN'%(i, j))
    
            if my_files:
                for name in my_files:
                    f = open(name, 'r')
                    
                    first_line = f.readline()
                    
                    if len(first_line) == 68:
                        
                        lines = f.readlines()
    
                        for num,line in enumerate(lines):
                            if len(line.split()) > 2:
                                if line.split()[2] == 'x_j(%)':
                                    table_start_value = num + 1
    
                            if len(line.split()) > 1:
                                if line.split()[1] == '[N_base]':
                                    n_base_value = int(line.split()[0])
    
                        lines = lines[table_start_value : table_start_value + n_base_value + 1]
                        stripped_lines = []
                    
                        for line in lines:
                            line = line.rstrip()
                            stripped_lines.append(line)
                    
                        data = ascii.read(stripped_lines[:-1])
                        t = Table(data) 
                    
                        my_xj = np.array(t['col2'])
                        my_age = np.array(t['col5'])

                        normalization_values[i, j] = np.sum(((my_xj / np.sum(my_xj)) * np.log10(my_age)))
                    
    return normalization_values


# In[6]:

def generate_normalized_mass_ages(dim_one, dim_two):
    """Generate ages normalized with respect to mass for all pertinent pixels in a 
    host galaxy image. 
    
    This function reads all of the output files and finds the normalization parameter Mini_j, which 
    is initial mass fraction, and the column of ages for this host galaxy. Age distributions are calculated
    by normalizing the ages with respect to Mini_j. This is an important output from STARLIGHT because 
    previous studies have indicated an age-distribution correlation between supernovae and host galaxies.
    """

    mass_normalization_values = np.empty((dim_one, dim_two))
    mass_normalization_values[:] = np.nan

    for i in range(dim_one):
        for j in range(dim_two):

            my_files = glob.glob('/Users/MCilento/STARLIGHT/KUG0210-078 files/KUG0210-078._%s_%s.C11.gm.CCM.BN'%(i, j))

            if my_files:
                for name in my_files:
                    f = open(name, 'r')
                    lines = f.readlines()

                    for num, line in enumerate(lines):
                        if len(line.split()) > 3:
                            if line.split()[3] == 'Mini_j(%)':
                                table_start_value = num + 1

                        if len(line.split()) > 1:
                            if line.split()[1] == '[N_base]':
                                n_base_value = int(line.split()[0])

                    lines = lines[table_start_value : table_start_value + n_base_value + 1]
                      
                    data = ascii.read(lines)
                    t = Table(data) 
                    
                    my_Mini_j = np.array(t['col3'])
                    my_age = np.array(t['col5'])
                    
                    mass_normalization_values[i, j] = np.sum(((my_Mini_j / np.sum(my_Mini_j)) * np.log10(my_age)))
                        
    return mass_normalization_values

# In[38]:

###This cell provides an example of generating the synthetic spectra for a single pixel
##It is not necessary to create a function to do this

fig = plt.figure(figsize=(15,10))

f = open('./manga-8158-1901 files/manga-8158-1901_11_14.C11.gm.CCM.BN', 'r')
f2 = open('./manga-8158-1901 files/manga-8158-1901_18_25.C11.gm.CCM.BN', 'r')

lines = f.readlines()
lines2 = f2.readlines()

for num,line in enumerate(lines):
    if len(line.split()) > 1:
        if line.split()[1] == '[fobs_norm': 
            f_obs = float(line.split()[0])

    if len(line.split()) > 1:
        if line.split()[1] == '[Nl_obs]':
            wave_len = int(line.split()[0])
            table_start = num + 1

            lines = lines[table_start : table_start + wave_len + 1]

            data = ascii.read(lines)
            t = Table(data) 

            wave = np.array(t['col1'])
            flux_syn = np.array(t['col3']) * f_obs
            flux_obs = np.array(t['col2']) * f_obs

for num,line in enumerate(lines2):
    if len(line.split()) > 1:
        if line.split()[1] == '[fobs_norm': 
            f_obs2 = float(line.split()[0])

    if len(line.split()) > 1:
        if line.split()[1] == '[Nl_obs]':
            wave_len2 = int(line.split()[0])
            table_start2 = num + 1

            lines2 = lines2[table_start2 : table_start2 + wave_len2 + 1]

            data2 = ascii.read(lines2)
            t2 = Table(data2) 

            wave2 = np.array(t2['col1'])
            flux_syn2 = np.array(t2['col3']) * f_obs2
            flux_obs2 = np.array(t2['col2']) * f_obs2
    
plt.plot(wave, flux_syn, color='blue', label='pixel [11, 14]')
plt.plot(wave2, flux_syn2, color='red', label='pixel [18, 25]')
plt.title('Synthetic Spectra for Resolved Pixels', fontsize=30, y = 1.02)
plt.ylabel(r'$Flux_\lambda \hspace{0.5}[normalized]$', fontsize=26)
plt.xlabel(r'$Wavelength \hspace{0.5}(\AA)$', fontsize=26)
plt.xticks(fontsize="17"); plt.yticks(fontsize="17")
plt.legend(fontsize="16")

plt.savefig('Poster2.png')

plt.show()


# In[27]:

def main():
    """Post-process output data files from STARLIGHT.
    
    This function reads in the data from STARLIGHT output files.  All other functions 
    are called within main(), and all data sets are returned so that after main() is 
    called one can make necessary images and plots for host galaxy analysis. 
    """
    survey_name = input('survey_name (2MASS or PMAS): ')
    PATH = input('path name: ')

    params = gen_run_parameters(survey_name, PATH)
    
    dim_one, dim_two = generate_data_shape(params, PATH)
    
    #av_values, vo_values, vd_values = generate_av_ov_vd_values(dim_one, dim_two)
    
    normalization_values = generate_normalized_light_ages(dim_one, dim_two)
    
    #mass_normalization_values = generate_normalized_mass_ages(dim_one, dim_two)
    
    return normalization_values


# In[28]:

normalization_values = main()


# In[8]:

##Extinction Parameter
#plt.imshow(av_values, cmap='jet', origin='lower')
#plt.colorbar()
#plt.show()

##Velocity
#plt.imshow(vo_values, vmin = 0, cmap='jet', origin='lower')
#plt.colorbar()
#plt.show()

##Velocity Distribution
#plt.imshow(vd_values, cmap='jet', origin='lower')
#plt.colorbar()
#plt.show()


# In[29]:

## Plot of age normalized w.r.t x_j (light fraction)

#fig = plt.figure(figsize=(15,10))

plt.imshow(normalization_values, cmap='seismic', origin='lower')
plt.title('Host Galaxy Stellar Age Distribution', fontsize=30, y = 1.02)
plt.xticks(fontsize="17"); plt.yticks(fontsize="17")
cbar = plt.colorbar()
cbar.set_label('Logarithmic Scale (Young to Old)', fontsize=26)

plt.savefig('Poster.png')

plt.show()


# In[ ]:



