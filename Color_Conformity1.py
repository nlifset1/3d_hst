print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import astropy.constants
from scipy import spatial
import math
from astropy.table import Table
from astropy.coordinates.sky_coordinate import SkyCoord
import random
import os
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#bring in the data#
print 'getting data'
data = ascii.read('color_tabled.dat')

#making lists for the plot#
lst_r = 10**np.linspace(1.2,3.6,13)

lst_q1 = []
lst_all1 = []
lst_sf1 = []
lst_q2 = []
lst_all2 = []
lst_sf2 = []
lst_q3 = []
lst_all3 = []
lst_sf3 = []
lst_q4 = []
lst_all4 = []
lst_sf4 = []

mask_q1 = ((data['z'] > 0.5) & (data['z'] <= 1.0) & (((data['vj'] < 0.92) & (data['uv'] > 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] > (0.88*data['vj'] +0.49)))))
mask_sf1 = ((data['z'] > 0.5) & (data['z'] <= 1.0) & (((data['vj'] < 0.92) & (data['uv'] < 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] < (0.88*data['vj'] +0.49))) | (data['vj']>1.5)))
mask_all1 = ((data['z'] > 0.5) & (data['z'] <= 1.0))
mask_q2 = ((data['z'] > 1.0) & (data['z'] <= 1.5) & (((data['vj'] < 0.92) & (data['uv'] > 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] > (0.88*data['vj'] +0.49)))))
mask_sf2 = ((data['z'] > 1.0) & (data['z'] <= 1.5) & (((data['vj'] < 0.92) & (data['uv'] < 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] < (0.88*data['vj'] +0.49))) | (data['vj']>1.5)))
mask_all2 = ((data['z'] > 1.0) & (data['z'] <= 1.5))
mask_q3 = ((data['z'] > 1.5) & (data['z'] <= 2.0) & (((data['vj'] < 0.92) & (data['uv'] > 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] > (0.88*data['vj'] +0.49)))))
mask_sf3 = ((data['z'] > 1.5) & (data['z'] <= 2.0) & (((data['vj'] < 0.92) & (data['uv'] < 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] < (0.88*data['vj'] +0.49))) | (data['vj']>1.5)))
mask_all3 = ((data['z'] > 1.5) & (data['z'] <= 2.0))
mask_q4 = ((data['z'] > 2.0) & (data['z'] <= 2.5) & (((data['vj'] < 0.92) & (data['uv'] > 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] > (0.88*data['vj'] +0.49)))))
mask_sf4 = ((data['z'] > 2.0) & (data['z'] <= 2.5) & (((data['vj'] < 0.92) & (data['uv'] < 1.3)) | ((data['vj'] > 0.8) & (data['vj'] < 1.6) & (data['uv'] < (0.88*data['vj'] +0.49))) | (data['vj']>1.5)))
mask_all4 = ((data['z'] > 2.0) & (data['z'] <= 2.5))

q1 = data[mask_q1]
sf1 = data[mask_sf1]
all1 = data[mask_all1]
q2 = data[mask_q2]
sf2 = data[mask_sf2]
all2 = data[mask_all2]
q3 = data[mask_q3]
sf3 = data[mask_sf3]
all3 = data[mask_all3]
q4 = data[mask_q4]
sf4 = data[mask_sf4]
all4 = data[mask_all4]

for i in range(len(lst_r)):
    radius_label = 'q%s' % (i+1)
    temp = 0
    temp_counts = 0
    for gal in q1:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_q1.append(temp)

    temp = 0
    temp_counts = 0
    for gal in q2:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_q2.append(temp)

    temp = 0
    temp_counts = 0
    for gal in q3:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_q3.append(temp)

    temp = 0
    temp_counts = 0
    for gal in q4:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_q4.append(temp)

    temp = 0
    temp_counts = 0
    for gal in sf1:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_sf1.append(temp)

    temp = 0
    temp_counts =0
    for gal in sf2:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_sf2.append(temp)

    temp = 0
    temp_counts = 0
    for gal in sf3:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_sf3.append(temp)

    temp = 0
    temp_counts = 0
    for gal in sf4:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_sf4.append(temp)

    temp = 0
    temp_counts = 0
    for gal in all1:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_all1.append(temp)

    temp = 0
    temp_counts = 0
    for gal in all2:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_all2.append(temp)

    temp = 0
    temp_counts = 0
    for gal in all3:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_all3.append(temp)

    temp = 0
    temp_counts = 0
    for gal in all4:
        if gal[radius_label] != 5:
            temp += gal[radius_label]
            temp_counts += 1
    temp = temp/float(temp_counts)
    lst_all4.append(temp)

#plotting#

fig,pylab.axes = pylab.subplots(1, 4, sharex = True, sharey = True)

a1 = pylab.axes[0]
a2 = pylab.axes[1]
a3 = pylab.axes[2]
a4 = pylab.axes[3]

a1.set_aspect('equal')
a2.set_aspect('equal')
a3.set_aspect('equal')
a4.set_aspect('equal')

a1.set_adjustable('box-forced')
a2.set_adjustable('box-forced')
a3.set_adjustable('box-forced')
a4.set_adjustable('box-forced')

a1.plot(lst_r,lst_q1, marker='o', color='r')
a1.plot(lst_r,lst_all1, marker='o', color='k')
a1.plot(lst_r,lst_sf1, marker='o', color='b')

a2.plot(lst_r,lst_q2, marker='o', color='r')
a2.plot(lst_r,lst_all2, marker='o', color='k')
a2.plot(lst_r,lst_sf2, marker='o', color='b')

a3.plot(lst_r,lst_q3, marker='o', color='r')
a3.plot(lst_r,lst_all3, marker='o', color='k')
a3.plot(lst_r,lst_sf3, marker='o', color='b')

a4.plot(lst_r,lst_q4, marker='o', color='r', label='Quiescent Centrals')
a4.plot(lst_r,lst_all4, marker='o', color='k', label='All Cnetrals')
a4.plot(lst_r,lst_sf4, marker='o', color='b', label='Star-Forming Centrals')

pylab.suptitle('Quiescent Satellite Percentage per Aperture Radius', fontsize=17)
fig.text(0.35, 0.01, "Aperture Radius (kpc)", fontsize=18)
fig.text(0.05, 0.9, "Quiescent Satellite Percentage (Q. sat./total sat.)", rotation = "vertical", fontsize=15)

a1.set_title('0.5 < z < 1.0', fontsize=10)
a2.set_title('1.0 < z < 1.5', fontsize=10)
a3.set_title('1.5 < z < 2.0', fontsize=10)
a4.set_title('2.0 < z < 1.5', fontsize=10)

fig.subplots_adjust(wspace=0.1)

a4.legend(loc=1)

a1.set_xscale('log')
a2.set_xscale('log')
a3.set_xscale('log')
a4.set_xscale('log')

pylab.ylim([0,1])
pylab.xlim([10,5000])

pylab.ion()
pylab.show()

