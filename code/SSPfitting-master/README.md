# Project Overview

This directory provides code for analyzing IFU data cubes with the STARLIGHT 
spectral fitting package. <File Name> contains complete descriptions of the
project background, implementation, methodology, results, and analysis. 

The two main subdirectories **STARLIGHT Pre-Process** and **STARLIGHT 
Post-Process** contain the code used in all stages of project as well as
further descriptions of the individual stages.


## 1. STARLIGHT Pre-Process

Prior to performing SSP synthesis with STARLIGHT, one must adhere to 
a specific series of data pre-processing steps outlined in the STARLIGHT 
manual. This subdirectory provides the code to carry out such steps and 
generate properly formatted input files, while the README provides an in-depth 
explanation of the steps themselves and how to run the code.


## 2. STARLIGHT Post-Process

After performing SSP synthesis with STARLIGHT, it is necessary to parse the
output files to generate interesting results such as mean ages, mean
metallicities, stellar mass estimates, star formation rates, star formation 
histories. This subdirectory contains multiple functions, each of which
performs a different post-process task with the output files. The README 
provides descriptions of each function and the generated results. 