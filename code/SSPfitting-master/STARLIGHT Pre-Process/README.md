# STARLIGHT Pre-Process

The purpose of the script *Pre_Process.py* is to carry out the 
pre-processing steps outlined in the STARLIGHT Manual. These steps are as follows:

> 1. Flux-Calibrate the spectrum
> 2. Correct flux for galactic extinction
> 3. Shift the spectrum into rest-frame
> 4. Resample flux, flux error, and flag values uniformly

This program currently runs with two different data models: MANGA and PMAS. 


## Running *Pre-Process.py*
### Command Line Arguments

In order to run this program, one must pull the script *Pre-Process.py* from the 
STARLIGHT Pre-Process directory. 

To run a fits file through this program, place the above script and 
the fits file in the same local directory. Use terminal to navigate to that 
directory. Once in that directory, there is a series of command line arguments 
that must be entered for the program to run. These are as follows:

> 1. Redshift (z)
> 2. Extinction Parameter (a)
> 3. Survey Name (n)
> 4. Filename (f)
> 5. Path (p)

The first three arguments are specific to the data that will be fit with 
STARLIGHT. Currently, Survey Name must be either MANGA or PMAS. The fourth 
argument is unique to the user's needs. This will act as the name of the 
directory created by the program that houses all of the input files for 
STARLIGHT. An example of a filename choice is shown in the **Examples** section. 
The last argument is simply the complete fits file name.

In the command line, these arguments must be typed in a specific format. This
is as follows: 

**$ python Pre_Process.py -z (value) -a (value) -n (name) -f (name) -p (path)**

If desired, one can simply copy and paste this format and replace each (value),
(name), and (path) with the correct information. 

### Examples

In the STARLIGHT Pre-Process directory, there is a text file called 
*Data Model Information.* This contains the command line arguments for two 
fits files that have previously been used in *Pre-Process.py.* Notice that the filename
choice is based on the fits file name, or path. It is usually convenient to 
use a filename that indicates the galaxy name.

Below is the command line used to run *Pre-Process.py* with the fits file named
manga-8158-1901-LINCUBE.fits: 

**$ python Pre_Process.py -z 0.038393 -a 0.2604 -n MANGA -f manga-8158-1901 -p manga-8158-1901-LINCUBE.fits**

The output of this run will be a new folder named "\<filename\> files" in the 
working directory, which contains the input spectra for STARLIGHT. In the same
working directory will be a file named "\<filename)\>\_grid_.in" which contains
heading information and the corresponding file names which contain all of the
input spectra. 


**Note:** *Pre-Process.py* only produces the inputs files necessary for 
STARLIGHT to run. STARLIGHT will not run automatically, and if one wishes to run 
STARLIGHT, then refer to the manual. This is also located in the STARLIGHT 
Pre-Process directory. 