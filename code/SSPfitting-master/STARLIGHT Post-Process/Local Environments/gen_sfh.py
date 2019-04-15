#!/usr/bin/env python
# coding: utf-8

"""This script generates Star Formation History per SN Type without
interpolating (Gaussian Processes) in the first stage. The final result will
be the average of all cumulative distributions of star formation rates / age
such that one can plot the average cumulative distribution per type
and compare."""

import numpy as np
import glob

from astropy.io import ascii
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
import matplotlib.pyplot as plt


__author__ = 'Meghan Cilento'


# Define temporary global variables
table_path = 'tab2_names_consistent.txt'
output_files_path = '/Users/MCilento/Github/STARLIGHT/local_sp files/*.BN'


def generate_sn_info(table_path):
    """Return a data table with information about the supernovae

    This function generates a table specifying SN names, types, and redshifts

    Arguments:
        table_path (str): the path indicating the data table location

    Returns:
        The data table specifying SN names, types, and redshifts
    """

    info_table = ascii.read(table_path, include_names=['SN', 'Type', 'z'])
    return info_table


info_table = generate_sn_info(table_path)


def generate_calculation_values(output_files_path, info_table):
    """Return the average star formation rates per SN type per age, an array
    of 23 ages, present solar mass values per SN environment, and processed
    solar masses per SN environment

    This function uses the STARLIGHT output files for 244 local SN
    environments to generate the information used to plot the star formation
    history per SN type

    Arguments:
        output_files_path (str) : the path indicating the STARLIGHT output
                                  files location
        info_table        (str) : the data table specifying SN names, types,
                                  and redshifts

    Returns:
        type_(II, IIn, Ic, IIb, Ia, Ibc, Ib)_values  (lists) :
        lists containing the average star formation rates per SN type for 23 ages
        sorted_ages                                  (array) :
        an array containing the 23 ages for which star formation rate is
        calculated

    """
    all_output_files = glob.glob(output_files_path)  # full paths of .BN files

    type_II_values = []
    type_IIn_values = []
    type_Ic_values = []
    type_IIb_values = []
    type_Ia_values = []
    type_Ibc_values = []
    type_Ib_values = []

    present_masses = []
    processed_masses = []

    for path in all_output_files:
        # Changing from absolute to relative file paths will make
        # splitting by '/' no longer necessary.
        supernova_name = path.split('/')[6].split('.')[0]

        with open(path) as ofile:
            lines = ofile.readlines()

            # Check to make sure that the output file contains data
            # If it does not, skip that file
            if len(lines) < 5:
                continue

            for num, line in enumerate(lines):
                line_list = line.split()
                if 'Mini_j(%)' in line_list:
                    table_start_value = num + 1

                elif '[N_base]' in line_list:
                    n_base_value = int(line_list[0])

                elif '[Mini_tot' in line_list:
                    mini_tot = float(line_list[0])

                elif '[Mcor_tot' in line_list:
                    mcor_tot = float(line_list[0])

            # Calculate the total mass of the SN local environment
            redshift = np.asarray(info_table['z'])[
                np.where(np.asarray(info_table['SN']) == supernova_name)]
            type = np.asarray(info_table['Type'])[
                np.where(np.asarray(info_table['SN']) == supernova_name)]
            solar_luminosity = 3.826e33  # Unit: erg / sec
            luminosity_distance = cosmo.luminosity_distance(
                redshift) * 3.0857e24  # Unit: cm
            # mini_tot unit: cm^2 * sec * solar mass / erg
            mass_local_env = mcor_tot * 4 * np.pi * (
                    luminosity_distance ** 2) / solar_luminosity  # Unit: Solar Mass
            processed_mass_local_env = mini_tot * 4 * np.pi * (
                    luminosity_distance ** 2) / solar_luminosity  # Unit: Solar Mass

            present_masses.append(mass_local_env.value)
            processed_masses.append(processed_mass_local_env.value)

            # Define arrays for mass percentages and ages
            lines = lines[
                    table_start_value: table_start_value + n_base_value + 1]
            data = ascii.read(lines)
            t = Table(data)
            t.keep_columns(['col3', 'col5'])
            percent_mass_t = (np.array(t['col3']))
            ages = np.array(t['col5'])

            # Remove repeated ages so there are only 23
            sorted_ages = []
            for age in ages:
                if age not in sorted_ages:
                    sorted_ages.append(age)
            sorted_ages.sort(reverse=True)
            sorted_ages = np.array(sorted_ages)

            # Define array for sorted mass percentages (total percent per age)
            sorted_percent_mass_t = []
            for a in sorted_ages:
                sorted_percent_mass_t.append(
                    np.sum(percent_mass_t[np.where(ages == a)]))
            sorted_percent_mass_t = np.array(sorted_percent_mass_t)

            # Calculate the age bin widths
            age_widths = []
            for i, age in enumerate(sorted_ages):
                if i == 0 and (i + 1) <= 22:
                    age_widths.append(((age + sorted_ages[i + 1]) / 2))
                if (i >= 1) and ((i + 1) <= 22):
                    age_widths.append(abs(((age + sorted_ages[i + 1]) / 2) - (
                        ((sorted_ages[i - 1] + age) / 2))))
            age_widths = np.array(age_widths)
            age_widths = np.append(age_widths, 2080000)

            ## Calculate the Star Formation Rate at all 23 ages
            SFR_t = (sorted_percent_mass_t / 100) * processed_mass_local_env / age_widths * 10000

            # Calculate the (non)normalized cumulative distribution of SFR
            cum_sum = np.cumsum(SFR_t)
            norm_cum_sum = cum_sum / np.sum(SFR_t).value

        # Split based on SN Type
        if (type == 'II') or (type == 'IIP') or (type == 'IIL'):
            #type_II_values.append(norm_cum_sum)
            type_II_values.append(cum_sum)

        if (type == 'Ia') or (type == 'Ia-pec') or (type == 'Ia-91bg') or (
                type == 'Ia-91T') or (type == 'Ia-02cx') or (type == 'I'):
            #type_Ia_values.append(norm_cum_sum)
            type_Ia_values.append(cum_sum)

        if (type == 'Ibc') or (type == 'Ibc-pec') or (type == 'Ic-BL') or (
                type == 'IIb') or (type == 'Ib') or (type == 'Ic'):
            #type_Ibc_values.append(norm_cum_sum)
            type_Ibc_values.append(cum_sum)

        if type == ('Ib'):
            #type_Ib_values.append(norm_cum_sum)
            type_Ib_values.append(cum_sum)

        if type == ('Ic'):
            #type_Ic_values.append(norm_cum_sum)
            type_Ic_values.append(cum_sum)

        if type == ('IIn'):
            #type_IIn_values.append(norm_cum_sum)
            type_IIn_values.append(cum_sum)

        if type == ('IIb'):
            #type_IIb_values.append(norm_cum_sum)
            type_IIb_values.append(cum_sum)

    return type_II_values, type_IIn_values, type_Ic_values, type_IIb_values, \
           type_Ia_values, type_Ibc_values, type_Ib_values, sorted_ages


type_II_values, type_IIn_values, type_Ic_values, type_IIb_values, \
type_Ia_values, type_Ibc_values, type_Ib_values, sorted_ages = \
    generate_calculation_values(output_files_path, info_table)



def plot_distributions(type_II_values, type_IIn_values, type_Ic_values,
                       type_IIb_values,
                       type_Ia_values, type_Ibc_values, type_Ib_values,
                       sorted_ages):

    type_II_averages = sum(type_II_values) / len(type_II_values)
    type_Ia_averages = sum(type_Ia_values) / len(type_Ia_values)
    type_Ibc_averages = sum(type_Ibc_values) / len(type_Ibc_values)
    type_IIn_averages = sum(type_IIn_values) / len(type_IIn_values)
    type_IIb_averages = sum(type_IIb_values) / len(type_IIb_values)
    type_Ic_averages = sum(type_Ic_values) / len(type_Ic_values)
    type_Ib_averages = sum(type_Ib_values) / len(type_Ib_values)

    #plt.plot(np.log10(sorted_ages), type_Ic_averages, color='g',
    #         label='Type Ic', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_Ib_averages, color='y',
    #         label='Type Ib', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_IIn_averages, color='m',
    #         label='Type IIn', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_II_averages, color='k',
    #         label='Type II', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_Ibc_averages, color='r',
    #         label='Type Ibc + IIb', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_IIb_averages, color='C1',
    #         label='Type IIb', linewidth=0.8)
    #plt.plot(np.log10(sorted_ages), type_Ia_averages, color='C0',
    #         label='Type Ia', linewidth=0.8)

    #Comparing star formation rates to the values for SNe II
    plt.plot(np.log10(sorted_ages), type_Ic_averages / type_II_averages, color='g', label='Type Ic', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_Ib_averages / type_II_averages, color='y', label='Type Ib', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_IIn_averages / type_II_averages, color='m', label='Type IIn', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_Ibc_averages / type_II_averages, color='r', label='Type Ibc + IIb', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_IIb_averages / type_II_averages, color='C1', label='Type IIb', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_Ia_averages / type_II_averages, color='C0', label='Type Ia', linewidth=0.8)
    plt.plot(np.log10(sorted_ages), type_II_averages / type_II_averages, color='k', label='Type II', linestyle='--', linewidth=1)

    plt.ylabel('Fraction of Mass')
    #plt.ylabel(r'Fraction of Mass [$SN_{X}$ / $SN_{II}$]')
    #plt.ylabel('Average SFR [$({M}_\odot\ yr^{-1}) * 10^{-4}$]', fontsize='11')
    plt.xlabel(r'log$_{10}$(age[yr])')
    #plt.axvline(x=7.31, color='k', linewidth=0.5, linestyle='--')
    #plt.axvline(x=7.61, color='k', linewidth=0.5, linestyle='--')
    #plt.axvline(x=8, color='k', linewidth=0.5, linestyle='--')
    plt.xlim(6, 10.25)
    #plt.ylim(0.7, 1.6)
    plt.legend(fontsize='9')
    plt.title('Star Formation History Per SN Type')
    plt.tight_layout()

    plt.savefig('SFH_II.png')

    plt.show()


plot_distributions(type_II_values, type_IIn_values, type_Ic_values,
                   type_IIb_values,
                   type_Ia_values, type_Ibc_values, type_Ib_values,
                   sorted_ages)




#Plot Star Formation Rate and Normalized Cumulative Distribution on the same set of axes.

#fig, ax1 = plt.subplots()
#
#color = 'tab:red'
#ax1.set_xlabel(r'log$_{10}$(age [yr])', fontsize='11')
#ax1.set_ylabel('Star Formation Rate [$({M}_\odot\ yr^{-1}) * 10^{-4}$]', color=color, fontsize='11')
#ax1.plot(np.log10(sorted_ages), SFR_t, color=color, label='Star Formation Rate')
#ax1.tick_params(axis='y', labelcolor=color)
#ax1.set_ylim([0, 176])
#
#ax2 = ax1.twinx()
#ax2.set_xlim([6, 10.25])
#
#color = 'tab:blue'
#ax2.set_ylabel('Fraction of Mass', color=color, fontsize='11')
#ax2.plot(np.log10(sorted_ages), norm_cum_dist, color=color, label='Normalized Cumulative Distribution')
#ax2.tick_params(axis='y', labelcolor=color)
#ax2.set_ylim([0, 1])
#
#plt.title('SN 1997cw: Type Ia-91T')
#fig.tight_layout()
#
#
#plt.savefig('SFR and Normalize Cumulative Distribution.png')
#
#plt.show()
