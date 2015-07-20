print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math
import random
from scipy import spatial
import os
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#making a median function easier#
def median(lst):
    return numpy.median(numpy.array(lst))

#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0))
data_flagged = data[idx]


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z_peak'] >= 0.5) & (gal['z_peak'] <= 2.5)):
            lst_gal_1.append([gal['id'], gal['field']])

#getting list of galaxies from previous list that also avoid edges by 0.05 degrees#
lst_gal =[]
upper_RA = 215.25
lower_RA = 34.27
upper_DEC = 62.33
lower_DEC = -28.0
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    if ((gal_info['ra'] < upper_RA) and (gal_info['ra'] > lower_RA)):
        if ((gal_info['dec'] < upper_DEC) and (gal_info['dec'] > lower_DEC)):
            lst_gal.append(gal)


#splitting the massed list into four redshift bins#
#bin 1: 2.0 < z < 2.5#
#bin 2: 1.5 < z < 2.0#
#bin 3: 1.0 < z < 1.5#
#bin 4: 0.5 < z < 2.0#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
for gal in lst_gal:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    if gal_info['z_peak'] >= 2.0:
        lst_gal1.append(gal)
    elif gal_info['z_peak'] >= 1.5:
        lst_gal2.append(gal)
    elif gal_info['z_peak'] >= 1.0:
        lst_gal3.append(gal)
    elif gal_info['z_peak'] >= 0.5:
        lst_gal4.append(gal)
        
lst_mass1 = []
lst_mass2 = []
lst_mass3 = []
lst_mass4 = []

for gal in lst_gal1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_mass1.append(gal_info['lmass'])
for gal in lst_gal2:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_mass2.append(gal_info['lmass'])
for gal in lst_gal3:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_mass3.append(gal_info['lmass'])
for gal in lst_gal4:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_mass4.append(gal_info['lmass'])

max1 = max(lst_mass1)
max2 = max(lst_mass2)
max3 = max(lst_mass3)
max4 = max(lst_mass4)

print 'bin1: %s, bin2: %s, bin3: %s, bin4: %s' % (max1[0], max2[0], max3[0], max4[0])



                            













