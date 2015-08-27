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
stuff = ascii.read('values_R2.dat')

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]

#initializing lists#
#breakdown: lst_errors [central mass number] s [satellite mass number]#
lst_errors1s1 = []
lst_errors2s1 = []
lst_errors3s1 = []
lst_errors4s1 = []
lst_errors1s2 = []
lst_errors2s2 = []
lst_errors3s2 = []
lst_errors4s2 = []
lst_errors1s3 = []
lst_errors2s3 = []
lst_errors3s3 = []
lst_errors4s3 = []
lst_errors1s4 = []
lst_errors2s4 = []
lst_errors3s4 = []
lst_errors4s4 = []

lst_radius = 10**np.linspace(1.2,3.6,13)

lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []

#binning based on mass#
print 'binning'
for gal in stuff:
    if gal['lmass'] <= 11.15:
        lst_gal1.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.3:
        lst_gal2.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.45:
        lst_gal3.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.8:
        lst_gal4.append([gal['id'],gal['field']])

#for each radius and each bin, calculating error bars via std of randoms#
for i in range(len(lst_radius)):
    print '%s' % (i)
    lst_rand1s1 = []
    lst_rand1s2 = []
    lst_rand1s3 = []
    lst_rand1s4 = []
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand1s1.append(gal_info['rand1r{}'.format(i+1)])
        lst_rand1s2.append(gal_info['rand2r{}'.format(i+1)])
        lst_rand1s3.append(gal_info['rand3r{}'.format(i+1)])
        lst_rand1s4.append(gal_info['rand4r{}'.format(i+1)])
    lst_errors1s1.append(np.std(lst_rand1s1))
    lst_errors1s2.append(np.std(lst_rand1s2))
    lst_errors1s3.append(np.std(lst_rand1s3))
    lst_errors1s4.append(np.std(lst_rand1s4))

    lst_rand2s1 = []
    lst_rand2s2 = []
    lst_rand2s3 = []
    lst_rand2s4 = []
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand2s1.append(gal_info['rand1r{}'.format(i+1)])
        lst_rand2s2.append(gal_info['rand2r{}'.format(i+1)])
        lst_rand2s3.append(gal_info['rand3r{}'.format(i+1)])
        lst_rand2s4.append(gal_info['rand4r{}'.format(i+1)])
    lst_errors2s1.append(np.std(lst_rand2s1))
    lst_errors2s2.append(np.std(lst_rand2s2))
    lst_errors2s3.append(np.std(lst_rand2s3))
    lst_errors2s4.append(np.std(lst_rand2s4))

    lst_rand3s1 = []
    lst_rand3s2 = []
    lst_rand3s3 = []
    lst_rand3s4 = []
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand3s1.append(gal_info['rand1r{}'.format(i+1)])
        lst_rand3s2.append(gal_info['rand2r{}'.format(i+1)])
        lst_rand3s3.append(gal_info['rand3r{}'.format(i+1)])
        lst_rand3s4.append(gal_info['rand4r{}'.format(i+1)])
    lst_errors3s1.append(np.std(lst_rand3s1))
    lst_errors3s2.append(np.std(lst_rand3s2))
    lst_errors3s3.append(np.std(lst_rand3s3))
    lst_errors3s4.append(np.std(lst_rand3s4))

    lst_rand4s1 = []
    lst_rand4s2 = []
    lst_rand4s3 = []
    lst_rand4s4 = []
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand4s1.append(gal_info['rand1r{}'.format(i+1)])
        lst_rand4s2.append(gal_info['rand2r{}'.format(i+1)])
        lst_rand4s3.append(gal_info['rand3r{}'.format(i+1)])
        lst_rand4s4.append(gal_info['rand4r{}'.format(i+1)])
    lst_errors4s1.append(np.std(lst_rand4s1))
    lst_errors4s2.append(np.std(lst_rand4s2))
    lst_errors4s3.append(np.std(lst_rand4s3))
    lst_errors4s4.append(np.std(lst_rand4s4))

#writing the table#
data = Table([lst_errors1s1,lst_errors1s2,lst_errors1s3,lst_errors1s4,lst_errors2s1,lst_errors2s2,lst_errors2s3,lst_errors2s4,lst_errors3s1,lst_errors3s2,lst_errors3s4,lst_errors3s4,lst_errors4s1,lst_errors4s2,lst_errors4s3,lst_errors4s4], names=['std1s1','std1s2','std1s3','std1s4','std2s1','std2s2','std2s3','std2s4','std3s1','std3s2','std3s3','std3s4','std4s1','std4s2','std4s3','std4s4'])
ascii.write(data, 'lmass_errors2.dat')



