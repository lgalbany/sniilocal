"""This script plots average light contribution vs. age for all
relevant supernova types in PISCO."""

import numpy as np
import glob
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

table = ascii.read("tab2.txt", include_names=['SN', 'Type'])

sn_names = np.asarray(table['SN'])
sn_types = np.asarray(table['Type'])


# Generate all output files
my_files = glob.glob('/Users/MCilento/Github/STARLIGHT/local_sp files/*.BN')

type_II_values = []
type_IIb_values = []
type_IIn_values = []
type_Ia_values = []
type_Ib_values = []
type_Ic_values = []
type_Ibc_values = []

for path in my_files:

    with open(path) as ofile:
        lines = ofile.readlines()

        if len(lines) < 5:
            continue

        for num, line in enumerate(lines):
            if len(line.split()) > 3:
                if line.split()[3] == 'Mini_j(%)':
                    table_start_value = num + 1

            if len(line.split()) > 1:
                if line.split()[1] == '[N_base]':
                    n_base_value = int(line.split()[0])

        lines = lines[table_start_value: table_start_value + n_base_value + 1]

        data = ascii.read(lines)
        t = Table(data)

        t.keep_columns(['col2', 'col5', 'col6'])

        # Define arrays for normalized light fraction, ages, and metallicity
        my_xj = (np.array(t['col2']))
        normalized_xj = (my_xj / np.sum(my_xj)) * 100

        my_ages = np.array(t['col5'])
        my_Z_j = np.array(t['col6'])

        # Define empty lists for sorted ages and light fraction totals per age
        sorted_ages = []
        x_j_totals = []

        for age in my_ages:
            if age not in sorted_ages:
                sorted_ages.append(age)

        sorted_ages.sort()

        final_ages = np.asarray(sorted_ages)

        for a in final_ages:
            x_j_totals.append(np.sum(normalized_xj[np.where(my_ages == a)]))

        x_j_totals = np.asarray(x_j_totals)

    path = path[48:].split('.')
    type = sn_types[np.where(sn_names == path[0])]

    new_ages = final_ages

    if type == 'II':
        type_II_values.append(x_j_totals)
    elif type == 'Ia':
        type_Ia_values.append(x_j_totals)
    elif type == 'Ib':
        type_Ib_values.append(x_j_totals)
    elif type == 'Ic':
        type_Ic_values.append(x_j_totals)
    elif type == 'Ibc':
        type_Ibc_values.append(x_j_totals)
    elif type == 'IIb':
        type_IIb_values.append(x_j_totals)
    elif type == 'IIn':
        type_IIn_values.append(x_j_totals)

plot_ages = np.log10(new_ages)

type_II_averages = sum(type_II_values) / len(type_II_values)
type_Ia_averages = sum(type_Ia_values) / len(type_Ia_values)
type_Ib_averages = sum(type_Ib_values) / len(type_Ib_values)
type_Ic_averages = sum(type_Ic_values) / len(type_Ic_values)
type_Ibc_averages = sum(type_Ibc_values) / len(type_Ibc_values)
type_IIb_averages = sum(type_IIb_values) / len(type_IIb_values)
type_IIn_averages = sum(type_IIn_values) / len(type_IIn_values)

plt.plot(plot_ages, type_II_averages, color='r', label='Type II')
plt.plot(plot_ages, type_Ia_averages, color='b', label='Type Ia')
plt.plot(plot_ages, type_Ib_averages, color='g', label='Type Ib')
plt.plot(plot_ages, type_Ic_averages, color='c', label='Type Ic')
plt.plot(plot_ages, type_Ibc_averages, color='k', label='Type Ibc')
plt.plot(plot_ages, type_IIb_averages, color='m', label='Type IIb')
plt.plot(plot_ages, type_IIn_averages, color='y', label='Type IIn')
plt.ylabel('Average Light Fraction Contribution [%]')
plt.xlabel('log10(age) [yr]')
plt.legend()
plt.show()
