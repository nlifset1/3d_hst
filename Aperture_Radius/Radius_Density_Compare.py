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

print 'bringing in data'
#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")
stuff = ascii.read('values_R.dat')

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z_peak'] >= 0.5) and (gal['z_peak'] <= 2.5)):
            lst_gal_1.append([gal['id'], gal['field']])

print 'edging'
#EDGING#
lst_gal = []
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    if gal[1] == 'AEGIS':
        if (gal_info['ra'] >= 3.746/(math.pi/180)) and (gal_info['ra'] <= 3.756821/(math.pi/180)):
            if (gal_info['dec'] >= 0.920312/(math.pi/180)) and (gal_info['dec'] <= 0.925897/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'COSMOS':
        if (gal_info['ra'] >= 2.619737/(math.pi/180)) and (gal_info['ra'] <= 2.620718/(math.pi/180)):
            if (gal_info['dec'] >= 0.038741/(math.pi/180)) and (gal_info['dec'] <= 0.043811/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'GOODS-N':
        if (gal_info['ra'] >= 3.298072/(math.pi/180)) and (gal_info['ra'] <= 3.307597/(math.pi/180)):
            if (gal_info['dec'] >= 1.084787/(math.pi/180)) and (gal_info['dec'] <= 1.087936/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'GOODS-S':
        if (gal_info['ra'] >= 0.925775/(math.pi/180)) and (gal_info['ra'] <= 0.929397/(math.pi/180)):
            if (gal_info['dec'] >= -0.487098/(math.pi/180)) and (gal_info['dec'] <= -0.483591/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'UDS':
        if (gal_info['ra'] >= 0.59815/(math.pi/180)) and (gal_info['ra'] <= 0.602889/(math.pi/180)):
            if (gal_info['dec'] >= -0.091376/(math.pi/180)) and (gal_info['dec'] <= -0.090305/(math.pi/180)):
                lst_gal.append(gal)

print 'listing'
#making lists for the plots of radius vs density, r is in KPC#
lst_r = 10**np.linspace(1.2,3.6,13)
lst_density = []
lst_rand = []
lst_final = []
lst_error = []
#for each radius, calculate average density with background subtraction#
for i in range(len(lst_r)):
    print '%s' % (i)
    radius_label = 'radius%s' % (i+1)
    rand_label = 'rand%s' % (i+1)
    density_total = 0
    density_rand_total = 0
    lst_special = []
    for gal in lst_gal:
        #finding number density of each galaxy at given radius and averaging them#
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        within = float(gal_info[radius_label])
        within_rand = float(gal_info[rand_label])
        lst_special.append(within_rand)
        density_total += within
        density_rand_total += within_rand
    density_ave = float(density_total)/len(lst_gal)
    density_rand_ave = float(density_rand_total)/len(lst_gal)
    final_ave = density_ave - density_rand_ave
    lst_density.append(density_ave)
    lst_rand.append(density_rand_ave)
    lst_final.append(final_ave)
    lst_error.append(np.std(lst_special))

#plotting radius vs density#
pylab.errorbar(lst_r, lst_final, yerr=lst_error, marker = 'o', markeredgecolor='none', color='black', linestyle='-', label='Selected Galaxies After Subraction')
pylab.plot(lst_r, lst_density, marker = 'o', markeredgecolor='none', color='b', linestyle='-', label = 'Selected Massive Galaxies')
pylab.errorbar(lst_r, lst_rand, marker = 'o', markeredgecolor='none', color='r', linestyle='-', label = 'Average Background Number Density')

pylab.suptitle('Galaxy Number Density per Aperture Radius at All Redshifts', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.0002])
pylab.yscale('log')
pylab.xscale('log')
pylab.legend(loc=3)


pylab.ion()
pylab.show()

