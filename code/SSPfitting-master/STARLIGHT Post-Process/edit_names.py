"""This script standardizes the names of galaxies and supernovae in PISCO.
There are multiple consistency checks to test the changes made to the names,
and they can be uncommented depending on what must be checked."""

from astropy.table import Table
from astropy.io import ascii
import glob
import numpy as np

"""##############################################################
Consistency Checks
##############################################################"""

# Make one big table containing a flag column indicating the 28 bad SN
original_table = Table.read('tab2_names_consistent.txt', format='ascii')
sn_names = np.asarray(original_table['SN'])

bad_sns = []

sns_removed = Table.read('tab2_28SNs_removed.txt', format='ascii')
cut_sn_names = np.array(sns_removed['SN'])

for name in sn_names:
    if name not in cut_sn_names:
        bad_sns.append(name)

print(bad_sns)

flag_list = []
for sn in sn_names:
    if sn in bad_sns:
        flag_list.append(1)
    elif sn not in bad_sns:
        flag_list.append(0)

flag = np.asarray(flag_list)

original_table['Flag'] = flag
new_table = original_table['GALname', 'GALID', 'Morphology', 'RAgalaxy',
                         'DECgalaxy', 'z', 'SN', 'Type', 'Flag', 'RASN',
                         'DECSN', 'Sep', 'Sample', 'Comments']

new_table['GALID'].description = "Galaxy ID in the CALIFA extended sample"
new_table['Flag'].description = "Flags indicating 28 SN with bad " \
                                "STARLIGHT fits. 1 = Bad. 0 = Good."
new_table['Sep'].description = "Supernova offset from galaxy core"
new_table['Sample'].description = ("Source. Sample meaning: \n"
                                  " 1- CALIFA DR3; \n"
                                  " 2- CALIFA extension projects (not included in DR3); \n"
                                  " 3- CALIFA dedicated project PI: Marino SNe: 2012ApJ...754...61M); \n"
                                  " 4- CALIFA pilot project (programme H09-3.5-068; 2012A%26A...545A..58S); \n"
                                  " 5- PPAK IFS Nearby Galaxies Survey (PINGS; 2010MNRAS.405..735R); \n"
                                  " 6- CALIFA programme H15-3.5-004; \n"
                                  " 7- CALIFA programme F16-3.5-006; \n"
                                  " 8- CALIFA programmes H16-3.5-012,F17-3.5-001, and H17-3.5-001. \n")
new_table.write("sn_targets.ecsv", format='ascii.ecsv', overwrite=True)


# Check the galaxies & SNe after removing the 28 bad SNe and changing naming
# convention.
#sns_removed = Table.read('tab2_28SNs_removed.txt', format='ascii')
#table_2 = Table.read('tab2_names_consistent.txt', format='ascii')
#
#all_files = glob.glob('/Users/MCilento/Github/STARLIGHT/total_sp/*.BN')
#
#directory_galaxy_names = []
#for path in all_files:
#    gal_name = path.split('/')[6].split('.')[0]
#    directory_galaxy_names.append(gal_name)
#
#for row in sns_removed:
#    row['GALname'] = row['GALname'].replace(' ', '').upper()
#    row['SN'] = row['SN'].replace(' ', '')
#
#sns_removed_new_names = Table.read('tab2_28SNs_removed_names_consistent.txt',
#                                   format='ascii')
#
#original_with_28 = list(table_2['GALname'])
#modified_no_28 = list(sns_removed_new_names['GALname'])
#
#for x in modified_no_28:
#    if x not in directory_galaxy_names:
#        print(x)


# Check all galaxies in local directory and modified galaxy names from tab2.txt
#all_files = glob.glob('/Users/MCilento/Github/STARLIGHT/total_sp/*.BN')
#
#directory_galaxy_names = []
#for path in all_files:
#    gal_name = path.split('/')[6].split('.')[0]
#    directory_galaxy_names.append(gal_name)
#
#
#table_2 = Table.read('tab2_names_consistent.txt', format='ascii')
#
#for row in table_2:
#    row['GALname'] = row['GALname'].replace(' ', '').upper()
#    row['SN'] = row['SN'].replace(' ', '')
#
#table_galaxy_names = list(table_2['GALname'])
#
#for x in table_galaxy_names:
#    if x not in directory_galaxy_names:
#        print(x)


"""##############################################################
Now we know that the local directory containing 232 galaxies is consistent
with the same 232 galaxies named in tab2.txt after the names have been 
modified with our convention:

Local directory names == tab2_names_consistent.txt

We also know that after removing the 28 SN and modifying their names with our 
convention the names in tab2_28SNs_removed_names_consistent.txt matches the 
names in the local directory:

Local directory names == tab2_28SNs_removed_names_consistent.txt

Going forward, it is best to use the two text files above that match the local
directory containing the spectrum files for these galaxies. Along with the
galaxy name modifications, the white spaces have been removed from the SN 
name, but nothing else was modified.
##############################################################"""


