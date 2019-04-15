#!/usr/bin/env python
# coding: utf-8

"""This code generates two different galaxy mass estimates from STARLIGHT
output files. STARLIGHT provides direct formulae for calculating both the
present galaxy mass in stars (solar mass units) as well as how many solar
masses have been processed into stars throughout the galaxy’s lifetime.
"""

import numpy as np
import glob

from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
from decimal import Decimal
from astropy.io import ascii


__author__ = 'Meghan Cilento'


# Define temporary global variables
# Later change these to relative file paths instead of absolute file paths
table_path = '/Users/MCilento/Github/SSPfitting/STARLIGHT Post-Process/tab2_names_consistent.txt'
output_files_path = '/Users/MCilento/Github/STARLIGHT/total_sp/*.BN'


def generate_table_info(table_path):
    """Return a data table with information about 232 PISCO galaxies

    Arguments:
        table_path (str): the path indicating the data table location

    Returns:
        The data table with PISCO galaxy information
    """

    info_table = Table.read(table_path, format='ascii')
    return info_table

info_table = generate_table_info(table_path)



def generate_calculation_values(output_files_path):
    """Return lists with galaxy names and the mass parameters needed to
    make mass estimates

    Arguments:
        output_files_path (str): the path indicating the STARLIGHT output
                                 files' location
    Returns:
        A list of all galaxy names' whose mass must be calculated

        A list of mini_tot values needed to calculate the processed solar
        masses over a galaxy's lifetime

        A list of mcor_tot values needed to calculate the present solar
        masses in a galaxy
     """

    # full paths of .BN files
    all_output_files = glob.glob(output_files_path)

    galaxy_names = []
    mini_tot_values = []
    mcor_tot_values = []

    for path in all_output_files:
        # Changing from absolute to relative file paths will make
        # splitting by '/' no longer necessary.
        galaxy_name = path.split('/')[6].split('.')[0]
        with open(path) as ofile:

            lines = ofile.readlines()

            # Check to make sure that the output file contains data
            # If it does not, skip that file
            if len(lines) < 5:
                continue

            galaxy_names.append(galaxy_name)

            for num, line in enumerate(lines):
                line_list = line.split()

                if '[Mini_tot' in line_list:
                    mini_tot = float(line_list[0])
                    mini_tot_values.append(mini_tot)

                if '[Mcor_tot' in line_list:
                    mcor_tot = float(line_list[0])
                    mcor_tot_values.append(mcor_tot)

    return galaxy_names, mini_tot_values, mcor_tot_values

galaxy_names, mini_tot_values, mcor_tot_values = generate_calculation_values(output_files_path)



def generate_mass_estimates(galaxy_names, info_table, mini_tot_values, mcor_tot_values):
    """Return a table containing the galaxy names, redshifts, present masses,
    and processed solar masses over a galaxy's lifetime

    Arguments:
        galaxy_names    (list) : A list of all galaxy names' whose mass must be
                                 calculated
        info_table      (table): The data table with PISCO galaxy information
        mini_tot_values (list) : A list of mini_tot values needed to calculate the processed solar
                                 masses over a galaxy's lifetime
        mcor_tot_values (list) : A list of mcor_tot values needed to calculate the present solar
                                 masses in a galaxy

    Returns:
        A table containing the galaxy names, redshifts, present masses, and
        processed solar masses over a galaxy's lifetime
        """

    # solar_luminosity units : erg / sec
    solar_luminosity = 3.826e33

    redshifts = []
    all_present_galaxy_mass = []
    all_processed_solar_masses = []

    for i, name in enumerate(galaxy_names):
        redshift = np.asarray(info_table['z'])[np.where(np.asarray(info_table['GALname']) == name)[0]]
        luminosity_distance = cosmo.luminosity_distance(redshift) * 3.0857e24  # Unit: cm

        # mini_tot_values units        : cm^2 * sec * solar mass / erg
        # mcor_tot_values units        : cm^2 * sec * solar mass / erg
        # present_galaxy_mass units    : solar mass
        # processed_solar_masses units : solar mass
        present_galaxy_mass = mcor_tot_values[i] * 4 * np.pi * (luminosity_distance ** 2) / solar_luminosity
        processed_solar_masses = mini_tot_values[i] * 4 * np.pi * (luminosity_distance ** 2) / solar_luminosity

        redshifts.append(redshift[0])
        all_present_galaxy_mass.append('%.3E' % Decimal(present_galaxy_mass[0].value))
        all_processed_solar_masses.append('%.3E' % Decimal(processed_solar_masses[0].value))

    table = Table({'Galaxy Name': galaxy_names,
                   'Redshift': redshifts,
                   'STARLIGHT Present Galaxy Mass': all_present_galaxy_mass,
                   'STARLIGHT Processed Solar Masses': all_processed_solar_masses},
                  names=['Galaxy Name', 'Redshift', 'Present Galaxy Mass',
                         'Processed Solar Masses'])

    return table

table = generate_mass_estimates(galaxy_names, info_table, mini_tot_values, mcor_tot_values)

# Write the table to a file
ascii.write(table, 'galaxy_masses.txt', overwrite=True, format='fixed_width')
