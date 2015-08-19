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
print 'getting data'
data = ascii.read('values_color.dat')

#getting error data#
errors = ascii.read('color_errors.dat')
lst_errorsrsr = []
lst_errorsrsb = []
lst_errorsbsr = []
lst_errorsbsb = []
for error in errors:
    lst_errorsrsr.append(error[0])
    lst_errorsrsb.append(error[1])
    lst_errorsbsr.append(error[2])
    lst_errorsbsb.append(error[3])

#splitting the central galaxies into two bins based on ssfr#
lst_galr = []
lst_galb = []
print 'binning'
for gal in data:
    uv = gal['uv']
    vj = gal['vj']
    if (((vj < 0.92) & (uv > 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv > (0.88*vj +0.49)))):
        lst_galr.append([gal['id'], gal['field']])
    elif (((vj < 0.92) & (uv < 1.3)) | ((vj > 0.8) & (vj < 1.6) & (uv < (0.88*vj +0.49))) | (vj>1.5)):
        lst_galb.append([gal['id'], gal['field']])


#making lists for the plot#
lst_r = 10**np.linspace(1.2,3.6,13)


lst_densityrsr = []
lst_densityrsb = []
lst_densitybsr = []
lst_densitybsb = []

for i in range(len(lst_r)):
    radius_labelr = 'densityrr%s' % (i+1)
    radius_labelb = 'densitybr%s' % (i+1)
    rand_labelr = 'randrr%s' % (i+1)
    rand_labelb = 'randbr%s' % (i+1)
    print '%s' % (radius_labelr)

    density_totalrsr = 0
    density_totalrsb = 0
    for gal in lst_galr:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalrsr += (densityr - density_randr)
        density_totalrsb += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_aversr = float(density_totalrsr)/len(lst_galr)
    density_aversb = float(density_totalrsb)/len(lst_galr)
    lst_densityrsr.append(density_aversr)
    lst_densityrsb.append(density_aversb)

    density_totalbsr = 0
    density_totalbsb = 0
    for gal in lst_galb:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalbsr += (densityr - density_randr)
        density_totalbsb += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_avebsr = float(density_totalbsr)/len(lst_galb)
    density_avebsb = float(density_totalbsb)/len(lst_galb)
    lst_densitybsr.append(density_avebsr)
    lst_densitybsb.append(density_avebsb)
    

    

#plotting radius vs density#

fig,pylab.axes = pylab.subplots(1, 2, sharex = True, sharey = True)

a1 = pylab.axes[0]
a2 = pylab.axes[1]

a1.errorbar(lst_r, lst_densityrsr, yerr=lst_errorsrsr, color='r', marker='o', markeredgecolor='none', linestyle='-', label = 'Quiescent Satellites')
a1.errorbar(lst_r, lst_densityrsb, yerr=lst_errorsrsb, color='b', marker='o', markeredgecolor='none', linestyle='-', label = 'Star-Forming Satellites')

a2.errorbar(lst_r, lst_densitybsr, yerr=lst_errorsbsr, color='r', marker='o', markeredgecolor='none', linestyle='-', label = 'Quiescent Satellites')
a2.errorbar(lst_r, lst_densitybsb, yerr=lst_errorsbsb, color='b', marker='o', markeredgecolor='none', linestyle='-', label = 'Star-Forming Satellites')

pylab.suptitle('Galaxy Number Density per Aperture Radius in Star-Forming Bins', fontsize=17)
fig.text(0.35, 0.01, "Aperture Radius (kpc)", fontsize=18)
fig.text(0.01, 0.85, "Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)", rotation = "vertical", fontsize=18)

a1.legend(loc=1, prop={'size':12})
a2.legend(loc=1, prop={'size':12})

a1.set_title('Quiescent Central Galaxies', fontsize=12)
a2.set_title('Star-Forming Central Galaxies', fontsize=12)

pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.002])

a1.set_yscale('log')
a2.set_yscale('log')

a1.set_xscale('log')
a2.set_xscale('log')

fig.subplots_adjust(hspace=0.1, wspace=0.1)

pylab.ion()
pylab.show()
