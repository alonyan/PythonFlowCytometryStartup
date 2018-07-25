
#%% load packages
import os
os.chdir('/Users/Alonyan/Documents/Python/example')

import numpy as np
import fit
import matplotlib.pyplot as plt
import glob
import flow_data as flow

# ======================================
# ======================================
## load data
# ======================================
# ======================================

#%% get data file names.  each txt file represents a unique dose of
# jak inhibitor
file_names = glob.glob('data/*.txt')

# load each text file into a unique flow class instances.
# classes, in a nutshell, allow a variable instance to
# have properties and a set of
# definied functions that operate on these properties

# instantiate a python list.  a list is a python data type
TestData = []
for wfile in file_names:
    # the '+=' syntax is similar to C, it appends the data list
    # with the new class instance
    TestData += [flow.data(wfile)]

#%% Call the print_channels function of the 0^th class instance
TestData[0].print_channels()
# the channel properties should print in command window, e.g. ipython terminal

## get pstat5 data from text files
pstat = np.array([])
for wdata in TestData:
    tmp_mean = np.mean(np.log10(wdata.get_data('pSTAT5')))
    pstat = np.hstack([pstat, tmp_mean])

## note an exponent is not a ^ in python, it is **
pstat = 10**pstat
inhib_dose = 1000 / 2**np.linspace(0, 11, 12)

#%% ======================================
# ======================================
## fit data to hill model
# ======================================
# ======================================

# fit simulated data
pars_mle = fit.fit_hill(inhib_dose, pstat)

# x data for plotting smooth model fit line
xplot = np.logspace(np.log10(inhib_dose.min()),
                    np.log10(inhib_dose.max()), 250)

#%% plot data and model
plt.figure()
plt.plot(inhib_dose, pstat, 'o',
         ms=10, mec='none',
         color='k',
         alpha=0.5)
plt.plot(xplot, fit.hill_model(pars_mle, xplot), '-',
         linewidth=2.5, color='k', alpha=0.7)
         
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'[I$_\mathregular{AZD1480}$] nM', fontsize=15)
plt.ylabel('pSTAT5 (a.u.)', fontsize=15)
plt.show(block=False)
