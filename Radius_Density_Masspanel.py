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
stuff = ascii.read('values_R2.dat')

#getting the lists for errors of each bin#
errors = ascii.read('lmass_errors2.dat')
lst_errors1s1 = []
lst_errors2s1 = []
lst_errors3s1 = []
lst_errors4s1 = []
lst_errors1s2 = []
lst_errors2s2 = []
lst_errors3s2 = []
lst_errors4s2 = []
lst_errors1s3 = []
lst_errors2s3 = []
lst_errors3s3 = []
lst_errors4s3 = []
lst_errors1s4 = []
lst_errors2s4 = []
lst_errors3s4 = []
lst_errors4s4 = []
for error in errors:
    lst_errors1s1.append(error[0])
    lst_errors1s2.append(error[1])
    lst_errors1s3.append(error[2])
    lst_errors1s4.append(error[3])
    lst_errors2s1.append(error[4])
    lst_errors2s2.append(error[5])
    lst_errors2s3.append(error[6])
    lst_errors2s4.append(error[7])
    lst_errors3s1.append(error[8])
    lst_errors3s2.append(error[9])
    lst_errors3s3.append(error[10])
    lst_errors3s4.append(error[11])
    lst_errors4s1.append(error[12])
    lst_errors4s2.append(error[13])
    lst_errors4s3.append(error[14])
    lst_errors4s4.append(error[15])

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


lst_density1s1 = []
lst_density2s1 = []
lst_density3s1 = []
lst_density4s1 = []
lst_density1s2 = []
lst_density2s2 = []
lst_density3s2 = []
lst_density4s2 = []
lst_density1s3 = []
lst_density2s3 = []
lst_density3s3 = []
lst_density4s3 = []
lst_density1s4 = []
lst_density2s4 = []
lst_density3s4 = []
lst_density4s4 = []


for i in range(len(lst_r)):
    
    radius_label = 'radius%s' % (i+1)
    rand_label = 'rand%s' % (i+1)
    print radius_label
    
    density_total1s1 = 0
    density_total1s2 = 0
    density_total1s3 = 0
    density_total1s4 = 0
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density1 = gal_info['density1r%s' % (i+1)]
        density2 = gal_info['density2r%s' % (i+1)]
        density3 = gal_info['density3r%s' % (i+1)]
        density4 = gal_info['density4r%s' % (i+1)]
        density_rand1 = gal_info['rand1r%s' % (i+1)]
        density_rand2 = gal_info['rand2r%s' % (i+1)]
        density_rand3 = gal_info['rand3r%s' % (i+1)]
        density_rand4 = gal_info['rand4r%s' % (i+1)]
        density_total1s1 += (density1 - density_rand1)
        density_total1s2 += (density2 - density_rand2)
        density_total1s3 += (density3 - density_rand3)
        density_total1s4 += (density4 - density_rand4)
    #averaging density of each galaxy at each radius#
    density_ave1s1 = float(density_total1s1)/len(lst_gal1)
    density_ave1s2 = float(density_total1s2)/len(lst_gal1)
    density_ave1s3 = float(density_total1s3)/len(lst_gal1)
    density_ave1s4 = float(density_total1s4)/len(lst_gal1)
    lst_density1s1.append(density_ave1s1)
    lst_density1s2.append(density_ave1s2)
    lst_density1s3.append(density_ave1s3)
    lst_density1s4.append(density_ave1s4)


    density_total4s1 = 0
    density_total4s2 = 0
    density_total4s3 = 0
    density_total4s4 = 0
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density1 = gal_info['density1r%s' % (i+1)]
        density2 = gal_info['density2r%s' % (i+1)]
        density3 = gal_info['density3r%s' % (i+1)]
        density4 = gal_info['density4r%s' % (i+1)]
        density_rand1 = gal_info['rand1r%s' % (i+1)]
        density_rand2 = gal_info['rand2r%s' % (i+1)]
        density_rand3 = gal_info['rand3r%s' % (i+1)]
        density_rand4 = gal_info['rand4r%s' % (i+1)]
        density_total1s1 += (density1 - density_rand1)
        density_total1s2 += (density2 - density_rand2)
        density_total1s3 += (density3 - density_rand3)
        density_total1s4 += (density4 - density_rand4)
    #averaging density of each galaxy at each radius#
    density_ave4s1 = float(density_total4s1)/len(lst_gal4)
    density_ave4s2 = float(density_total4s2)/len(lst_gal4)
    density_ave4s3 = float(density_total4s3)/len(lst_gal4)
    density_ave4s4 = float(density_total4s4)/len(lst_gal4)
    lst_density4s1.append(density_ave4s1)
    lst_density4s2.append(density_ave4s2)
    lst_density4s3.append(density_ave4s3)
    lst_density4s4.append(density_ave4s4)

    density_total3s1 = 0
    density_total3s2 = 0
    density_total3s3 = 0
    density_total3s4 = 0
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density1 = gal_info['density1r%s' % (i+1)]
        density2 = gal_info['density2r%s' % (i+1)]
        density3 = gal_info['density3r%s' % (i+1)]
        density4 = gal_info['density4r%s' % (i+1)]
        density_rand1 = gal_info['rand1r%s' % (i+1)]
        density_rand2 = gal_info['rand2r%s' % (i+1)]
        density_rand3 = gal_info['rand3r%s' % (i+1)]
        density_rand4 = gal_info['rand4r%s' % (i+1)]
        density_total1s1 += (density1 - density_rand1)
        density_total1s2 += (density2 - density_rand2)
        density_total1s3 += (density3 - density_rand3)
        density_total1s4 += (density4 - density_rand4)
    #averaging density of each galaxy at each radius#
    density_ave3s1 = float(density_total3s1)/len(lst_gal3)
    density_ave3s2 = float(density_total3s2)/len(lst_gal3)
    density_ave3s3 = float(density_total3s3)/len(lst_gal3)
    density_ave3s4 = float(density_total3s4)/len(lst_gal3)
    lst_density3s1.append(density_ave3s1)
    lst_density3s2.append(density_ave3s2)
    lst_density3s3.append(density_ave3s3)
    lst_density3s4.append(density_ave3s4)


    density_total2s1 = 0
    density_total2s2 = 0
    density_total2s3 = 0
    density_total2s4 = 0
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        density1 = gal_info['density1r%s' % (i+1)]
        density2 = gal_info['density2r%s' % (i+1)]
        density3 = gal_info['density3r%s' % (i+1)]
        density4 = gal_info['density4r%s' % (i+1)]
        density_rand1 = gal_info['rand1r%s' % (i+1)]
        density_rand2 = gal_info['rand2r%s' % (i+1)]
        density_rand3 = gal_info['rand3r%s' % (i+1)]
        density_rand4 = gal_info['rand4r%s' % (i+1)]
        density_total1s1 += (density1 - density_rand1)
        density_total1s2 += (density2 - density_rand2)
        density_total1s3 += (density3 - density_rand3)
        density_total1s4 += (density4 - density_rand4)
    #averaging density of each galaxy at each radius#
    density_ave2s1 = float(density_total2s1)/len(lst_gal2)
    density_ave2s2 = float(density_total2s2)/len(lst_gal2)
    density_ave2s3 = float(density_total2s3)/len(lst_gal2)
    density_ave2s4 = float(density_total2s4)/len(lst_gal2)
    lst_density2s1.append(density_ave2s1)
    lst_density2s2.append(density_ave2s2)
    lst_density2s3.append(density_ave2s3)
    lst_density2s4.append(density_ave2s4)

    
#plotting radius vs density#

fig,pylab.axes = pylab.subplots(2, 2, sharex = True, sharey = True)

a1 = pylab.axes[0,0]
a2 = pylab.axes[0,1]
a3 = pylab.axes[1,0]
a4 = pylab.axes[1,1]

a1.errorbar(lst_r, lst_density1s1, yerr=lst_errors1s1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '9.415 < Sat. Lmass < 9.8')
a1.errorbar(lst_r, lst_density1s2, yerr=lst_errors1s2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '9.8 < Sat. LMass < 10.3')
a1.errorbar(lst_r, lst_density1s3, yerr=lst_errors1s3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '10.3 < Sat. LMass < 10.8')
a1.errorbar(lst_r, lst_density1s4, yerr=lst_errors1s4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '10.8 < Sat. LMass')

a2.errorbar(lst_r, lst_density2s1, yerr=lst_errors2s1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '9.415 < Sat. Lmass < 9.8')
a2.errorbar(lst_r, lst_density2s2, yerr=lst_errors2s2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '9.8 < Sat. LMass < 10.3')
a2.errorbar(lst_r, lst_density2s3, yerr=lst_errors2s3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '10.3 < Sat. LMass < 10.8')
a2.errorbar(lst_r, lst_density2s4, yerr=lst_errors2s4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '10.8 < Sat. LMass')

a3.errorbar(lst_r, lst_density3s1, yerr=lst_errors3s1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '9.415 < Sat. Lmass < 9.8')
a3.errorbar(lst_r, lst_density3s2, yerr=lst_errors3s2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '9.8 < Sat. LMass < 10.3')
a3.errorbar(lst_r, lst_density3s3, yerr=lst_errors3s3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '10.3 < Sat. LMass < 10.8')
a3.errorbar(lst_r, lst_density3s4, yerr=lst_errors3s4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '10.8 < Sat. LMass')

a4.errorbar(lst_r, lst_density4s1, yerr=lst_errors4s1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '9.415 < Sat. Lmass < 9.8')
a4.errorbar(lst_r, lst_density4s2, yerr=lst_errors4s2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '9.8 < Sat. LMass < 10.3')
a4.errorbar(lst_r, lst_density4s3, yerr=lst_errors4s3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '10.3 < Sat. LMass < 10.8')
a4.errorbar(lst_r, lst_density4s4, yerr=lst_errors4s4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '10.8 < Sat. LMass')


pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Mass Bins', fontsize=17)
fig.text(0.44, 0.01, "Aperture Radius (kpc)", fontsize=18)
fig.text(0.01, 0.9, "Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)", rotation = "vertical", fontsize=18)

a1.legend(loc=1, prop={'size':9})
a2.legend(loc=1, prop={'size':9})
a3.legend(loc=1, prop={'size':9})
a4.legend(loc=1, prop={'size':9})

pylab.xlim([10,5000])
pylab.ylim([0.0000001,0.002])

a1.set_yscale('log')
a2.set_yscale('log')
a3.set_yscale('log')
a4.set_yscale('log')

a1.set_xscale('log')
a2.set_xscale('log')
a3.set_xscale('log')
a4.set_xscale('log')

fig.subplots_adjust(hspace=0.1, wspace=0.1)


pylab.ion()
pylab.show()

