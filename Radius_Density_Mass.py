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
errors = ascii.read('errors.dat')
lst_errors = []
for error in errors:
    lst_errors.append(error)

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z_peak'] >= 0.5) and (gal['z_peak'] <= 2.5)):
            lst_gal_1.append([gal['id'], gal['field']])

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
    


#splitting the massed list into four lmass bins#
#bin 1: 11.0 < lmass < 11.15#
#bin 2: 11.15 < lmass < 11.3#
#bin 3: 11.3 < lmass < 11.45#
#bin 4: 11.45 < lmass < 11.6#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
for gal in lst_gal:
    gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
    if gal_info['lmass'] <= 11.15:
        lst_gal1.append(gal)
    elif gal_info['lmass'] <= 11.3:
        lst_gal2.append(gal)
    elif gal_info['lmass'] <= 11.45:
        lst_gal3.append(gal)
    elif gal_info['lmass'] <= 11.6:
        lst_gal4.append(gal)

#making lists for the plots of radius vs density, r is in KPC#
lst_r = [20,30,50,75,100,200,300,500,750,1000]
lst_density1 = []
lst_density2 = []
lst_density3 = []
lst_density4 = []
for r in lst_r:
    radius_label = 'r%s' % (r)

    density_total1 = 0
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_total1 += density
    #averaging density of each galaxy at each radius#
    density_ave1 = float(density_total1)/len(lst_gal1)
    lst_density1.append(density_ave1)

    density_total2 = 0
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_total2 += density
    #averaging density of each galaxy at each radius#
    density_ave2 = float(density_total2)/len(lst_gal2)
    lst_density2.append(density_ave2)

    density_total3 = 0
    for gal in lst_gal_massed3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_total3 += density
    #averaging density of each galaxy at each radius#
    density_ave3 = float(density_total3)/len(lst_gal3)
    lst_density3.append(density_ave3)

    density_total4 = 0
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_total4 += density
    #averaging density of each galaxy at each radius#
    density_ave4 = float(density_total4)/len(lst_gal4)
    lst_density4.append(density_ave4)


#plotting radius vs density#

pylab.plot(lst_r, lst_density1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '11.0 < LMass < 11.15')
pylab.plot(lst_r, lst_density2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '11.15 < LMass < 11.3')
pylab.plot(lst_r, lst_density3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '11.3 < LMass < 11.45')
pylab.plot(lst_r, lst_density4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '11.45 < LMass < 11.6')

pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Mass Bins', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.legend(loc=1)
pylab.xlim([15,1050])
pylab.ylim([0.0000005,0.002])
pylab.yscale('log')
pylab.xscale('log')



pylab.ion()
pylab.show()

