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

#splitting the massed list into four redshift bins#
#bin 1: 2.0 < z < 2.5#
#bin 2: 1.5 < z < 2.0#
#bin 3: 1.0 < z < 1.5#
#bin 4: 0.5 < z < 2.0#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
print 'binning'
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
lst_r = 10**np.linspace(1.2,3.6,13)

lst_density1 = []
lst_density2 = []
lst_density3 = []
lst_density4 = []

for i in range(len(lst_r)):
    radius_label = 'radius%s' % (i+1)
    rand_label = 'rand%s' % (i+1)
    print radius_label
    density_total1 = 0
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total1 += (density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave1 = float(density_total1)/len(lst_gal1)
    lst_density1.append(density_ave1)

    density_total2 = 0
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total2 += (density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave2 = float(density_total2)/len(lst_gal2)
    lst_density2.append(density_ave2)

    density_total3 = 0
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total3 += (density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave3 = float(density_total3)/len(lst_gal3)
    lst_density3.append(density_ave3)

    density_total4 = 0
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density = gal_info[radius_label]
        density_rand = gal_info[rand_label]
        density_total4 += (density - density_rand)
    #averaging density of each galaxy at each radius#
    density_ave4 = float(density_total4)/len(lst_gal4)
    lst_density4.append(density_ave4)

#making lists for Tomer's data to plot#
#1: 0.04<z<0.07#
#2: 0.4<z<0.9#
#3: 0.9<z<1.3#
#4: 1.3<z<1.6#
lst_tomer1 = [0.0,0.0,0.0,0.0,10**(-5.2),10**(-5.34),10**(-5.6),10**(-5.82),10**(-6.06),10**(-6.26),10**(-6.52),10**(-6.7),10**(-7.09)]
lst_tomer2 = [0.0,10**(-4.35),10**(-4.26),10**(-5.01),10**(-4.78),10**(-5.03),10**(-5.29),10**(-5.44),10**(-5.71),10**(-5.94),10**(-6.02),10**(-6.21),10**(-6.4)]
lst_tomer3 = [10**(-4.14),10**(-4.35),10**(-4.7),10**(-4.89),10**(-4.93),10**(-5.28),10**(-5.73),10**(-5.77),10**(-5.86),10**(-5.98),10**(-6.29),10**(-6.53),10**(-7.15)]
lst_tomer4 = [10**(-4.1),10**(-4.31),10**(-4.79),10**(-4.78),10**(-5.07),10**(-5.38),10**(-5.45),10**(-5.69),10**(-5.65),10**(-5.78),10**(-5.95),10**(-6.11),10**(-6.38)]

lst_terror1 = [0.0,0.0,0.0,0.0,10**(-5.69),10**(-6.30),10**(-6.32),10**(-6.52),10**(-6.63),10**(-7.09),10**(-6.99),10**(-7.17),10**(-7.37)]
lst_terror2 = [0.0,10**(-5.23),10**(-5.47),10**(-5.61),10**(-5.76),10**(-5.98),10**(-6.07),10**(-6.24),10**(-6.37),10**(-6.45),10**(-6.55),10**(-6.71),10**(-6.81)]
lst_terror3 = [10**(-5.16),10**(-5.36),10**(-5.6),10**(-5.71),10**(-5.97),10**(-6.11),10**(-6.31),10**(-6.42),10**(-6.56),10**(-6.71),10**(-6.88),10**(-6.88),10**(-7.01)]
lst_terror4 = [10**(-5.02),10**(-5.28),10**(-5.41),10**(-5.63),10**(-5.8),10**(-6.0),10**(-6.23),10**(-6.34),10**(-6.44),10**(-6.55),10**(-6.61),10**(-6.73),10**(-6.79)]

lst_tr = [10**(1.2),10**(1.4),10**(1.6),10**(1.8),10**(2.0),10**(2.2),10**(2.4),10**(2.6),10**(2.8),10**(3.0),10**(3.2),10**(3.4),10**(3.6)]

        
#plotting radius vs density#

pylab.errorbar(lst_r, lst_density1, yerr=lst_errors1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '2.0 < z < 2.5')
pylab.errorbar(lst_r, lst_density2, yerr=lst_errors2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '1.5 < z < 2.0')
pylab.errorbar(lst_r, lst_density3, yerr=lst_errors3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '1.0 < z < 1.5')
pylab.errorbar(lst_r, lst_density4, yerr=lst_errors4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '0.5 < z < 1.0')
pylab.errorbar(lst_tr, lst_tomer1, yerr=lst_terror1, color='c', marker='o', markeredgecolor='none', linestyle='-', label = 'Tal\'s 0.04 < z < 0.07')
pylab.errorbar(lst_tr, lst_tomer2, yerr=lst_terror2, color='m', marker='o', markeredgecolor='none', linestyle='-', label = 'Tal\'s 0.4 < z < 0.9')
pylab.errorbar(lst_tr, lst_tomer3, yerr=lst_terror3, color='k', marker='o', markeredgecolor='none', linestyle='-', label = 'Tal\'s 0.9 < z < 1.3')
pylab.errorbar(lst_tr, lst_tomer4, yerr=lst_terror4, color='0.75', marker='o', markeredgecolor='none', linestyle='-', label = 'Tal\'s 1.3 < z < 1.6')

pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Redshift Bins', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.legend(loc=1)
pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.002])
pylab.yscale('log')
pylab.xscale('log')



pylab.ion()
pylab.show()

