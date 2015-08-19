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
print 'getting data'
data = ascii.read('values_color.dat')


#CENTRALS#
print 'counting centrals'
gal_number = len(data)

mask_r = (((data['vj'] < 0.92) & (data['uv'] > 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] > (0.88*data['vj'] +0.49))))
mask_b = (((data['vj'] < 0.92) & (data['uv'] < 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] < (0.88*data['vj'] +0.49))) | (data['vj']>1.5))

data_r = data[mask_r]
data_b = data[mask_b]

R_gal_number = len(data_r)
B_gal_number = len(data_b)



#SATELLITES#
print 'counting satellites'
data_me = ascii.read("C:\\3d_hst/3dhst_big.dat")
data_3dhst = ascii.read("C:\\3d_hst/3dhst_master.phot.v4.1.cat")

mask = ((data_me['z'] <= 2.5) & (data_me['z'] >= 0.5) & (data_3dhst['lmass'] >= 9.415) & (data_3dhst['use_phot'] == 1))
data_flagged = data_me[mask]
data_other = data_3dhst[mask]

sat_number = len(data_flagged)
R_sat_number = 0
B_sat_number = 0

for i in range(len(data_flagged)):
    gal = data_flagged[i]
    uv = -2.5*np.log10(gal['L153']/gal['L155'])
    vj = -2.5*np.log10(gal['L155']/gal['L161'])
    if (((vj < 0.92) & (uv > 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv > (0.88*vj +0.49)))):
        R_sat_number += 1
    elif (((vj < 0.92) & (uv < 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv < (0.88*vj +0.49))) | (vj>1.5)):
        B_sat_number += 1

print 'tabling'

data = [[gal_number, R_gal_number, B_gal_number, sat_number, R_sat_number, B_sat_number]]
columns = ('Central Galaxies', 'Quiescent Centrals', 'Star-Forming Centrals', 'Potential Satellites', 'Quiescent Satellites', 'Star-Forming Satellites')

rows = ('number')

# Add a table at the bottom of the axes
the_table = plt.table(cellText=data, colLabels=columns, loc='center')
plt.axis('off')
the_table.set_fontsize(22)
the_table.scale(1.2, 2)

