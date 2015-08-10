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

#making a median function easier#
def median(lst):
    return np.median(np.array(lst))

def mad(lst):
    return median(np.absolute(median(lst)-lst))

#bring in the data#
print 'getting data'
data = ascii.read("3dhst_master.phot.v4.1.cat")
stuff = ascii.read('values_r.dat')

#getting the lists for errors of each bin#
errors = ascii.read('lmass_errors.dat')
lst_errors1 = []
lst_errors2 = []
lst_errors3 = []
lst_errors4 = []
for error in errors:
    lst_errors1.append(error[0])
    lst_errors2.append(error[1])
    lst_errors3.append(error[2])
    lst_errors4.append(error[3])

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]    


#splitting the massed list into four lmass bins#
#bin 1: 11.0 < lmass < 11.15#
#bin 2: 11.15 < lmass < 11.3#
#bin 3: 11.3 < lmass < 11.45#
#bin 4: 11.45 < lmass < 11.6#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
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

#making lists for the plots of radius vs density, r is in KPC#
lst_r = 10**np.linspace(1.2,3.6,13)


lst_density1 = []
lst_density2 = []
lst_density3 = []
lst_density4 = []

lst_sigma1 = []
lst_sigma2 = []
lst_sigma3 = []
lst_sigma4 = []

for i in range(len(lst_r)):
    
    radius_label = 'radius%s' % (i+1)
    rand_label = 'rand%s' % (i+1)
    print radius_label
    
    density_total1 = 0
    lst_dens1 = []
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total1 += (density - density_rand)
        lst_dens1.append(density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave1 = float(density_total1)/len(lst_gal1)
    lst_density1.append(density_ave1)
    lst_sigma1.append(mad(lst_dens1))

    density_total2 = 0
    lst_dens2 = []
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total2 += (density - density_rand)
        lst_dens2.append(density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave2 = float(density_total2)/len(lst_gal2)
    lst_density2.append(density_ave2)
    lst_sigma2.append(mad(lst_dens2))

    density_total3 = 0
    lst_dens3 = []
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total3 += (density - density_rand)
        lst_dens3.append(density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave3 = float(density_total3)/len(lst_gal3)
    lst_density3.append(density_ave3)
    lst_sigma3.append(mad(lst_dens3))

    density_total4 = 0
    lst_dens4 = []
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total4 += (density - density_rand)
        lst_dens4.append(density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave4 = float(density_total4)/len(lst_gal4)
    lst_density4.append(density_ave4)
    lst_sigma4.append(mad(lst_dens4))



#plotting radius vs density#

pylab.errorbar(lst_r, lst_density1, yerr=lst_errors1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '11.0 < LMass < 11.15')
pylab.errorbar(lst_r, lst_density2, yerr=lst_errors2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '11.15 < LMass < 11.3')
pylab.errorbar(lst_r, lst_density3, yerr=lst_errors3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '11.3 < LMass < 11.45')
pylab.errorbar(lst_r, lst_density4, yerr=lst_errors4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '11.45 < LMass')

pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Mass Bins', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.legend(loc=1)
pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.002])
pylab.yscale('log')
pylab.xscale('log')



pylab.ion()
pylab.show()

