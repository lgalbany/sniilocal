"""This script performs gaussian process regression to the plots of light
fraction vs. age for a single supernova local environment."""


import numpy as np
import glob
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.table import Table

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern


my_files = glob.glob('/Users/MCilento/Github/STARLIGHT/local_sp files/1997cw.txt.C11.gm.CCM.BN')

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

        # Define arrays for normalized light fractions and  ages
        my_xj = (np.array(t['col2']))
        normalized_xj = (my_xj / np.sum(my_xj)) * 100

        my_ages = np.array(t['col5'])

        # Define empty lists for sorted ages and light fraction totals per age
        sorted_ages = []
        xj_totals = []

        for age in my_ages:
            if age not in sorted_ages:
                sorted_ages.append(age)

        sorted_ages.sort()

        final_ages = np.asarray(sorted_ages)

        for a in final_ages:
            xj_totals.append(np.sum(normalized_xj[np.where(my_ages == a)]))


        """Apply Gaussian Process Regression
        x is the training data: the sampled x values
        y is the target values: the actual data points
        X is the known age values: the actual 23 data points
        """

        # x Values
        x = np.atleast_2d(np.linspace(6, 10.255, 1000)).T

        # X Values
        X = np.log10(final_ages)
        X = np.reshape(X, (23, 1))

        # y Values
        y = np.asarray(xj_totals)
        y = np.reshape(y, 23)


        # Instantiate a Gaussian Process model
        kernel = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=1.5)
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, y)

        # Make the prediction on the meshed x-axis
        y_pred, sigma = gp.predict(x, return_std=True)


        """New Tasks
        --Get rid of any negative values: default to zero
        --Do the integral of the fit and normalize
            (array of 1000 values / integral * 100)
        """

        print(np.shape(y_pred))

        plt.plot(X, y, 'r.', markersize=10, label = u'Observed Data')
        plt.plot(X, y, 'g--', linewidth=2)
        plt.plot(x, y_pred, 'b-', label = u'Fit: Matern')
        plt.xlabel('Age [yr]')
        plt.ylabel('Normalized Light Fraction [%]')
        plt.legend(loc='upper left')
        plt.show()


if __name__ == '__main__':
    print('hello world')
