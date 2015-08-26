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

lst_c1 = []
lst_c2 = []
lst_c3 = []
lst_c4 = []
for i in range(len(lst_r)):
    print i
    lst_c1.append((lst_q1[i] - lst_sf1[i])/lst_all1[i])
    lst_c2.append((lst_q2[i] - lst_sf2[i])/lst_all2[i])
    lst_c3.append((lst_q3[i] - lst_sf3[i])/lst_all3[i])
    lst_c4.append((lst_q4[i] - lst_sf4[i])/lst_all4[i])
    

#plotting#

pylab.plot(lst_r, lst_c1, color='r', marker='o', label='0.5 < z < 1.0')
pylab.plot(lst_r, lst_c2, color='g', marker='o', label='1.0 < z < 1.5')
pylab.plot(lst_r, lst_c3, color='b', marker='o', label='1.5 < z < 2.0')
pylab.plot(lst_r, lst_c4, color='y', marker='o', label='2.0 < z < 2.5')

pylab.suptitle('Conformity per Aperture Radius in Four Redshift Bins', fontsize=17)
pylab.xlabel("Aperture Radius (kpc)", fontsize=18)
pylab.ylabel("Conformity", fontsize=15)

pylab.legend(loc=1)

pylab.xscale('log')

pylab.ylim([0,4])
pylab.xlim([10,5000])

pylab.ion()
pylab.show()

