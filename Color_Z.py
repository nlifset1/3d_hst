from astropy.io import ascii
import numpy as np
import math
import random
import astropy.constants
from scipy import spatial
import os
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
from operator import add
import matplotlib.pyplot as plt
import pylab
os.chdir('C:\\3d_hst')

#bring in the data#
print 'Reading in catalog...'
data_me = ascii.read("C:\\3d_hst/3dhst_big.dat")
data = ascii.read("C:\\3d_hst/3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
print 'masking'
mask = ((data_me['z'] <= 2.5) & (data_me['z'] >= 0.5) & (data['lmass'] >= 11.0) & (data['use_phot'] == 1))
data_flagged = data_me[mask]
data_other = data[mask]
print 'there are %s galaxies' % (len(data_flagged))


#initializing lists and counts for quiescent and star-forming#
lst_r_uv1 = []
lst_r_vj1 = []
lst_b_uv1 = []
lst_b_vj1 = []
lst_r_uv2 = []
lst_r_vj2 = []
lst_b_uv2 = []
lst_b_vj2 = []
lst_r_uv3 = []
lst_r_vj3 = []
lst_b_uv3 = []
lst_b_vj3 = []
lst_r_uv4 = []
lst_r_vj4 = []
lst_b_uv4 = []
lst_b_vj4 = []
quiescent = 0
star_forming = 0

#calculating uv and vj for each galaxy and sorting the galaxies based on that info#
#r is for quiescent, b is for star-forming#
#also binning based on redshift, with low number for low redshift#
for i in range(len(data_flagged)):
    gal = data_flagged[i]
    uv = -2.5*np.log10(gal['L153']/gal['L155'])
    vj = -2.5*np.log10(gal['L155']/gal['L161'])
    if (((vj < 0.92) & (uv > 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv > (0.88*vj +0.49)))):
        if data_other[i]['z_peak'] < 1.0:
            lst_r_uv1.append(uv)
            lst_r_vj1.append(vj)
        elif data_other[i]['z_peak'] < 1.5:
            lst_r_uv2.append(uv)
            lst_r_vj2.append(vj)
        elif data_other[i]['z_peak'] < 2.0:
            lst_r_uv3.append(uv)
            lst_r_vj3.append(vj)
        elif data_other[i]['z_peak'] < 2.5:
            lst_r_uv4.append(uv)
            lst_r_vj4.append(vj)
        quiescent += 1
    elif (((vj < 0.92) & (uv < 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv < (0.88*vj +0.49))) | (vj>1.5)):
        if data_other[i]['z_peak'] < 1.0:
            lst_b_uv1.append(uv)
            lst_b_vj1.append(vj)
        elif data_other[i]['z_peak'] < 1.5:
            lst_b_uv2.append(uv)
            lst_b_vj2.append(vj)
        elif data_other[i]['z_peak'] < 2.0:
            lst_b_uv3.append(uv)
            lst_b_vj3.append(vj)
        elif data_other[i]['z_peak'] < 2.5:
            lst_b_uv4.append(uv)
            lst_b_vj4.append(vj)
        star_forming += 1


#plotting#
        
pylab.plot(lst_r_vj1, lst_r_uv1, color='#fcae91', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='0.5<z<1.0')
pylab.plot(lst_r_vj2, lst_r_uv2, color='#fb6a4a', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='1.0<z<1.5')
pylab.plot(lst_r_vj3, lst_r_uv3, color='#de2d26', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='1.5<z<2.0')
pylab.plot(lst_r_vj4, lst_r_uv4, color='#a50f15', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='2.0<z<2.5')

pylab.plot(lst_b_vj1, lst_b_uv1, color='#9ecae1', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='0.5<z<1.0')
pylab.plot(lst_b_vj2, lst_b_uv2, color='#6baed6', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='1.0<z<1.5')
pylab.plot(lst_b_vj3, lst_b_uv3, color='#3182bd', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='1.5<z<2.0')
pylab.plot(lst_b_vj4, lst_b_uv4, color='#08519c', markeredgecolor='none', marker='o', linestyle='none', alpha=0.9, label='2.0<z<2.5')

pylab.text(0.05,2.35, 'Quiescent(%s)' % (quiescent), fontsize=12, fontweight='bold')
pylab.text(0.5,0.08, 'Star-Forming(%s)' % (star_forming), fontsize=12, fontweight='bold')

pylab.plot([0, 0.92], [1.3, 1.3], 'k-')
pylab.plot([0.92, 1.6], [1.3, 1.898], 'k-')
pylab.plot([1.6, 1.6], [1.898, 2.5], 'k-')


        
pylab.xlim([0,2.5])
pylab.ylim([0,2.5])

pylab.xlabel('(V-J)$_{rest}$', fontsize=18)
pylab.ylabel('(U-V)$_{rest}$', fontsize=18)
pylab.suptitle('Quiescent and Star-Forming Central Galaxies in Four Redshift Bins', fontsize=18)

pylab.ion()
pylab.show()



