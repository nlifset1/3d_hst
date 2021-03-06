print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import math
import random
import os
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")
stuff = ascii.read('values_R.dat')

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]

#initializing lists#
#low number for HIGH redshift bin#
lst_errors1 = []
lst_errors2 = []
lst_errors3 = []
lst_errors4 = []

lst_radius = 10**np.linspace(1.2,3.6,13)

lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []

#binning the galaxies based on redshift#
print 'binning'
for gal in stuff:
    if gal['z'] >= 2.0:
        lst_gal1.append([gal['id'],gal['field']])
    elif gal['z'] >= 1.5:
        lst_gal2.append([gal['id'],gal['field']])
    elif gal['z'] >= 1.0:
        lst_gal3.append([gal['id'],gal['field']])
    elif gal['z'] >= 0.5:
        lst_gal4.append([gal['id'],gal['field']])

#for each radius and each bin, calculating error bars via std of randoms#
for i in range(len(lst_radius)):
    print '%s' % (i)
    lst_rand1 = []
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand1.append(gal_info['rand{}'.format(i+1)])
    lst_errors1.append(np.std(lst_rand1))

    lst_rand2 = []
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand2.append(gal_info['rand{}'.format(i+1)])
    lst_errors2.append(np.std(lst_rand2))

    lst_rand3 = []
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand3.append(gal_info['rand{}'.format(i+1)])
    lst_errors3.append(np.std(lst_rand3))

    lst_rand4 = []
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand4.append(gal_info['rand{}'.format(i+1)])
    lst_errors4.append(np.std(lst_rand4))

#writing the table#
data = Table([lst_errors1,lst_errors2,lst_errors3,lst_errors4], names=['std1','std2','std3','std4'])
ascii.write(data, 'z_errors.dat')



