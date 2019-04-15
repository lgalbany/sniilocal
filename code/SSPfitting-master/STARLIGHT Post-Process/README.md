# STARLIGHT Post-Process

This subdirectory contains multiple functions, each of which
performs a different post-process task with STARLIGHT output files. Depending 
on the desired application, it may be necessary to use all functions or simply 
one function. 

The functions have been separated into two different folders, **Galaxies**
and **Local Environments.** Below are general descriptions of the folders
and the functions contained within them. 


## Galaxies
This folder contains functions which produce global galaxy properties such as
mass estimates and stellar age distribution. 

#### *gen_gal_masses.py*
Used to generate two different galaxy mass estimates from STARLIGHT
output files. STARLIGHT provides direct formulae for calculating both the
present galaxy mass in stars (solar mass units) as well as how many solar
masses have been processed into stars throughout the galaxy’s lifetime.

#### *gen_age_dist.py*
Note: This needs to be updated. Currently this is used for general
post-processing purposes. Ultimately, it will be used to produce stellar age
distributions.


## Local Environments
This folder contains functions which produce local supernova environment
properties such as average light fraction contributions and star formation
history.

#### *l_frac_fit.py*
Used to perform gaussian process regression and generate plots of light
fraction contribution vs. age for a single supernova local environment.

#### *l_frac_sn_type.py*
Used to generate plots of average light contribution vs. age for all
relevant supernova types in PISCO.

#### *gen_sfh.py*
Used to generate star formation history per SN Type without
interpolating (Gaussian Processes) in the first stage. The final result is
be the average of all cumulative distributions of star formation rates per age
such that one can plot the average cumulative distribution per type
and compare to Type II SN.