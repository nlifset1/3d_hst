print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import math
import random
import os
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")
stuff = ascii.read('values_r.dat')

#getting the lists for errors of each bin#
errors = ascii.read('lmass_errors.dat')
lst_errors1 = []
lst_errors2 = []
lst_errors3 = []
lst_errors4 = []
for error in errors:
    lst_errors1.append(error[0][0])
    lst_errors2.append(error[1][0])
    lst_errors3.append(error[2][0])
    lst_errors4.append(error[3][0])

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]

#splitting the massed list into four redshift bins#
#bin 1: 2.0 < z < 2.5#
#bin 2: 1.5 < z < 2.0#
#bin 3: 1.0 < z < 1.5#
#bin 4: 0.5 < z < 2.0#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
for gal in stuff:
    if gal['z'] >= 2.0:
        lst_gal1.append([gal['id'],gal['field']])
    elif gal['z'] >= 1.5:
        lst_gal2.append([gal['id'],gal['field']])
    elif gal['z'] >= 1.0:
        lst_gal3.append([gal['id'],gal['field']])
    elif gal['z'] >= 0.5:
        lst_gal4.append([gal['id'],gal['field']])

#making lists for the plots of radius vs density, r is in KPC#
lst_r = [20,30,50,75,100,200,300,500,750,1000]

lst_density1 = []
lst_density2 = []
lst_density3 = []
lst_density4 = []

lst_dens1 = []
lst_dens2 = []
lst_dens3 = []
lst_dens4 = []

for r in lst_r:
    radius_label = 'r%s' % (r)
    density_total1 = 0
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        lst_dens1.append(density)
        density_total1 += density
    #averaging density of each galaxy at each radius#
    density_ave1 = float(density_total1)/len(lst_gal1)
    lst_density1.append(density_ave1)

    density_total2 = 0
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        lst_dens2.append(density)
        density_total2 += density
    #averaging density of each galaxy at each radius#
    density_ave2 = float(density_total2)/len(lst_gal2)
    lst_density2.append(density_ave2)

    density_total3 = 0
    for gal in lst_gal:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        lst_dens3.append(density)
        density_total3 += density
    #averaging density of each galaxy at each radius#
    density_ave3 = float(density_total3)/len(lst_gal3)
    lst_density3.append(density_ave3)

    density_total4 = 0
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        lst_dens4.append(density)
        density_total4 += density
    #averaging density of each galaxy at each radius#
    density_ave4 = float(density_total4)/len(lst_gal4)
    lst_density4.append(density_ave4)

lst_r1 = []
lst_r2 = []
lst_r3 = []
lst_r4 = []
for radius in lst_r:
    for gal in lst_gal1:
        lst_r1.append(radius)
for radius in lst_r:
    for gal in lst_gal2:
        lst_r2.append(radius)
for radius in lst_r:
    for gal in lst_gal3:
        lst_r3.append(radius)
for radius in lst_r:
    for gal in lst_gal4:
        lst_r4.append(radius)
        
#plotting radius vs density#

pylab.errorbar(lst_r, lst_density1, yerr=lst_errors1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '2.0 < z < 2.5')
pylab.errorbar(lst_r, lst_density2, yerr=lst_errors2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '1.5 < z < 2.0')
pylab.errorbar(lst_r, lst_density3, yerr=lst_errors3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '1.0 < z < 1.5')
pylab.errorbar(lst_r, lst_density4, yerr=lst_errors4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '0.5 < z < 1.0')

pylab.scatter(lst_r1, lst_dens1, color='r', alpha=0.3)
pylab.scatter(lst_r2, lst_dens2, color='b', alpha=0.3)
pylab.scatter(lst_r3, lst_dens3, color='g', alpha=0.3)
pylab.scatter(lst_r4, lst_dens4, color='y', alpha=0.3)

pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Redshift Bins', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.legend(loc=1)
pylab.xlim([15,1050])
pylab.ylim([0.0000005,0.002])
pylab.yscale('log')
pylab.xscale('log')



pylab.ion()
pylab.show()

