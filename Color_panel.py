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
lst_errorsrsr1 = []
lst_errorsrsb1 = []
lst_errorsbsr1 = []
lst_errorsbsb1 = []
lst_errorsrsr2 = []
lst_errorsrsb2 = []
lst_errorsbsr2 = []
lst_errorsbsb2 = []
lst_errorsrsr3 = []
lst_errorsrsb3 = []
lst_errorsbsr3 = []
lst_errorsbsb3 = []
lst_errorsrsr4 = []
lst_errorsrsb4 = []
lst_errorsbsr4 = []
lst_errorsbsb4 = []
for error in errors:
    lst_errorsrsr1.append(error[0])
    lst_errorsrsb1.append(error[1])
    lst_errorsbsr1.append(error[2])
    lst_errorsbsb1.append(error[3])
    lst_errorsrsr2.append(error[4])
    lst_errorsrsb2.append(error[5])
    lst_errorsbsr2.append(error[6])
    lst_errorsbsb2.append(error[7])
    lst_errorsrsr3.append(error[8])
    lst_errorsrsb3.append(error[9])
    lst_errorsbsr3.append(error[10])
    lst_errorsbsb3.append(error[11])
    lst_errorsrsr4.append(error[12])
    lst_errorsrsb4.append(error[13])
    lst_errorsbsr4.append(error[14])
    lst_errorsbsb4.append(error[15])

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

#splitting the central galaxies into more bins based on redshift#
#low number is low redshift#
lst_galr1 = []
lst_galb1 = []
lst_galr2 = []
lst_galb2 = []
lst_galr3 = []
lst_galb3 = []
lst_galr4 = []
lst_galb4 = []
for gal in lst_galr:
    gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
    if gal_info['z'] < 1.0:
        lst_galr1.append(gal)
    elif gal_info['z'] < 1.5:
        lst_galr2.append(gal)
    elif gal_info['z'] < 2.0:
        lst_galr3.append(gal)
    elif gal_info['z'] < 2.5:
        lst_galr4.append(gal)

for gal in lst_galb:
    gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
    if gal_info['z'] < 1.0:
        lst_galb1.append(gal)
    elif gal_info['z'] < 1.5:
        lst_galb2.append(gal)
    elif gal_info['z'] < 2.0:
        lst_galb3.append(gal)
    elif gal_info['z'] < 2.5:
        lst_galb4.append(gal)


#making lists for the plot#
lst_r = 10**np.linspace(1.2,3.6,13)



lst_densityrsr1 = []
lst_densityrsb1 = []
lst_densitybsr1 = []
lst_densitybsb1 = []
lst_densityrsr2 = []
lst_densityrsb2 = []
lst_densitybsr2 = []
lst_densitybsb2 = []
lst_densityrsr3 = []
lst_densityrsb3 = []
lst_densitybsr3 = []
lst_densitybsb3 = []
lst_densityrsr4 = []
lst_densityrsb4 = []
lst_densitybsr4 = []
lst_densitybsb4 = []

for i in range(len(lst_r)):
    radius_labelr = 'densityrr%s' % (i+1)
    radius_labelb = 'densitybr%s' % (i+1)
    rand_labelr = 'randrr%s' % (i+1)
    rand_labelb = 'randbr%s' % (i+1)
    print '%s' % (radius_labelr)

    density_totalrsr1 = 0
    density_totalrsb1 = 0
    for gal in lst_galr1:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalrsr1 += (densityr - density_randr)
        density_totalrsb1 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_aversr1 = float(density_totalrsr1)/len(lst_galr1)
    density_aversb1 = float(density_totalrsb1)/len(lst_galr1)
    lst_densityrsr1.append(density_aversr1)
    lst_densityrsb1.append(density_aversb1)

    density_totalbsr1 = 0
    density_totalbsb1 = 0
    for gal in lst_galb1:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalbsr1 += (densityr - density_randr)
        density_totalbsb1 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_avebsr1 = float(density_totalbsr1)/len(lst_galb1)
    density_avebsb1 = float(density_totalbsb1)/len(lst_galb1)
    lst_densitybsr1.append(density_avebsr1)
    lst_densitybsb1.append(density_avebsb1)

    density_totalrsr2 = 0
    density_totalrsb2 = 0
    for gal in lst_galr2:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalrsr2 += (densityr - density_randr)
        density_totalrsb2 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_aversr2 = float(density_totalrsr2)/len(lst_galr2)
    density_aversb2 = float(density_totalrsb2)/len(lst_galr2)
    lst_densityrsr2.append(density_aversr2)
    lst_densityrsb2.append(density_aversb2)

    density_totalbsr2 = 0
    density_totalbsb2 = 0
    for gal in lst_galb2:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalbsr2 += (densityr - density_randr)
        density_totalbsb2 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_avebsr2 = float(density_totalbsr2)/len(lst_galb2)
    density_avebsb2 = float(density_totalbsb2)/len(lst_galb2)
    lst_densitybsr2.append(density_avebsr2)
    lst_densitybsb2.append(density_avebsb2)

    density_totalrsr3 = 0
    density_totalrsb3 = 0
    for gal in lst_galr3:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalrsr3 += (densityr - density_randr)
        density_totalrsb3 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_aversr3 = float(density_totalrsr3)/len(lst_galr3)
    density_aversb3 = float(density_totalrsb3)/len(lst_galr3)
    lst_densityrsr3.append(density_aversr3)
    lst_densityrsb3.append(density_aversb3)

    density_totalbsr3 = 0
    density_totalbsb3 = 0
    for gal in lst_galb3:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalbsr3 += (densityr - density_randr)
        density_totalbsb3 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_avebsr3 = float(density_totalbsr3)/len(lst_galb3)
    density_avebsb3 = float(density_totalbsb3)/len(lst_galb3)
    lst_densitybsr3.append(density_avebsr3)
    lst_densitybsb3.append(density_avebsb3)

    density_totalrsr4 = 0
    density_totalrsb4 = 0
    for gal in lst_galr4:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalrsr4 += (densityr - density_randr)
        density_totalrsb4 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_aversr4 = float(density_totalrsr4)/len(lst_galr4)
    density_aversb4 = float(density_totalrsb4)/len(lst_galr4)
    lst_densityrsr4.append(density_aversr4)
    lst_densityrsb4.append(density_aversb4)

    density_totalbsr4 = 0
    density_totalbsb4 = 0
    for gal in lst_galb4:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        densityr = gal_info[radius_labelr]
        densityb = gal_info[radius_labelb]
        density_randr = gal_info[rand_labelr]
        density_randb = gal_info[rand_labelb]
        density_totalbsr4 += (densityr - density_randr)
        density_totalbsb4 += (densityb - density_randb)
    #averaging density of each galaxy at each radius#
    density_avebsr4 = float(density_totalbsr4)/len(lst_galb4)
    density_avebsb4 = float(density_totalbsb4)/len(lst_galb4)
    lst_densitybsr4.append(density_avebsr4)
    lst_densitybsb4.append(density_avebsb4)
    

    
#removing data points from list that equal 0#
lst_densitybsr1 = np.delete(lst_densitybsr1, [0,1])
lst_errorsbsr1 = np.delete(lst_errorsbsr1, [0,1])
lst_densityrsb2 = np.delete(lst_densityrsb2, 0)
lst_errorsrsb2 = np.delete(lst_errorsrsb2, 0)
lst_densitybsr2 = np.delete(lst_densitybsr2, [0,1,2])
lst_errorsbsr2 = np.delete(lst_errorsbsr2, [0,1,2])

#making special radius lists for density lists just made#
lst_r1 = np.delete(lst_r, [0,1])
lst_r2 = np.delete(lst_r, 0)
lst_r3 = np.delete(lst_r, [0,1,2])


#plotting radius vs density#

fig,pylab.axes = pylab.subplots(2, 4, sharex = True, sharey = True)

a1 = pylab.axes[0,0]
a2 = pylab.axes[0,1]
a3 = pylab.axes[0,2]
a4 = pylab.axes[0,3]
a5 = pylab.axes[1,0]
a6 = pylab.axes[1,1]
a7 = pylab.axes[1,2]
a8 = pylab.axes[1,3]

a1.errorbar(lst_r, lst_densityrsr1, yerr=lst_errorsrsr1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = 'Quiescent Satellites')
a1.errorbar(lst_r, lst_densityrsb1, yerr=lst_errorsrsb1, color='b', marker='o', markeredgecolor='none', linestyle='-', label = 'Star-Forming Satellites')

a5.errorbar(lst_r1, lst_densitybsr1, yerr=lst_errorsbsr1, color='r', marker='o', markeredgecolor='none', linestyle='-')
a5.errorbar(lst_r, lst_densitybsb1, yerr=lst_errorsbsb1, color='b', marker='o', markeredgecolor='none', linestyle='-')

a2.errorbar(lst_r, lst_densityrsr2, yerr=lst_errorsrsr2, color='r', marker='o', markeredgecolor='none', linestyle='-')
a2.errorbar(lst_r2, lst_densityrsb2, yerr=lst_errorsrsb2, color='b', marker='o', markeredgecolor='none', linestyle='-')

a6.errorbar(lst_r3, lst_densitybsr2, yerr=lst_errorsbsr2, color='r', marker='o', markeredgecolor='none', linestyle='-')
a6.errorbar(lst_r, lst_densitybsb2, yerr=lst_errorsbsb2, color='b', marker='o', markeredgecolor='none', linestyle='-')

a3.errorbar(lst_r, lst_densityrsr3, yerr=lst_errorsrsr3, color='r', marker='o', markeredgecolor='none', linestyle='-')
a3.errorbar(lst_r, lst_densityrsb3, yerr=lst_errorsrsb3, color='b', marker='o', markeredgecolor='none', linestyle='-')

a7.errorbar(lst_r, lst_densitybsr3, yerr=lst_errorsbsr3, color='r', marker='o', markeredgecolor='none', linestyle='-')
a7.errorbar(lst_r, lst_densitybsb3, yerr=lst_errorsbsb3, color='b', marker='o', markeredgecolor='none', linestyle='-')

a4.errorbar(lst_r, lst_densityrsr4, yerr=lst_errorsrsr4, color='r', marker='o', markeredgecolor='none', linestyle='-')
a4.errorbar(lst_r, lst_densityrsb4, yerr=lst_errorsrsb4, color='b', marker='o', markeredgecolor='none', linestyle='-')

a8.errorbar(lst_r, lst_densitybsr4, yerr=lst_errorsbsr4, color='r', marker='o', markeredgecolor='none', linestyle='-', label = 'Quiescent Satellites')
a8.errorbar(lst_r, lst_densitybsb4, yerr=lst_errorsbsb4, color='b', marker='o', markeredgecolor='none', linestyle='-', label = 'Star-Forming Satellites')

a1.set_aspect('equal', adjustable='box-forced')
a2.set_aspect('equal', adjustable='box-forced')
a3.set_aspect('equal', adjustable='box-forced')
a4.set_aspect('equal', adjustable='box-forced')
a5.set_aspect('equal', adjustable='box-forced')
a6.set_aspect('equal', adjustable='box-forced')
a7.set_aspect('equal', adjustable='box-forced')
a8.set_aspect('equal', adjustable='box-forced')



pylab.suptitle('Galaxy Number Density per Aperture Radius in Star-Forming Bins', fontsize=17)
fig.text(0.35, 0.01, "Aperture Radius (kpc)", fontsize=18)
fig.text(0.0, 0.85, "Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)", rotation = "vertical", fontsize=18)

a1.set_title('0.5 < z < 1.0', fontsize=10)
a2.set_title('1.0 < z < 1.5', fontsize=10)
a3.set_title('1.5 < z < 2.0', fontsize=10)
a4.set_title('2.0 < z < 1.5', fontsize=10)

a1.set_ylabel('Quiescent Centrals')
a5.set_ylabel('Star-Forming Centrals')

pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.0002])

a1.set_yscale('log')
a2.set_yscale('log')
a3.set_yscale('log')
a4.set_yscale('log')
a5.set_yscale('log')
a6.set_yscale('log')
a7.set_yscale('log')
a8.set_yscale('log')

a1.set_xscale('log')
a2.set_xscale('log')
a3.set_xscale('log')
a4.set_xscale('log')
a5.set_xscale('log')
a6.set_xscale('log')
a7.set_xscale('log')
a8.set_xscale('log')


fig.subplots_adjust(hspace=0.0001, wspace=0.0001)

pylab.ion()
pylab.show()
