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
#bin 4: 11.45 < lmass < 11.8#
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
        density_total4s1 += (density1 - density_rand1)
        density_total4s2 += (density2 - density_rand2)
        density_total4s3 += (density3 - density_rand3)
        density_total4s4 += (density4 - density_rand4)
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
        density_total3s1 += (density1 - density_rand1)
        density_total3s2 += (density2 - density_rand2)
        density_total3s3 += (density3 - density_rand3)
        density_total3s4 += (density4 - density_rand4)
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
        density_total2s1 += (density1 - density_rand1)
        density_total2s2 += (density2 - density_rand2)
        density_total2s3 += (density3 - density_rand3)
        density_total2s4 += (density4 - density_rand4)
    #averaging density of each galaxy at each radius#
    density_ave2s1 = float(density_total2s1)/len(lst_gal2)
    density_ave2s2 = float(density_total2s2)/len(lst_gal2)
    density_ave2s3 = float(density_total2s3)/len(lst_gal2)
    density_ave2s4 = float(density_total2s4)/len(lst_gal2)
    lst_density2s1.append(density_ave2s1)
    lst_density2s2.append(density_ave2s2)
    lst_density2s3.append(density_ave2s3)
    lst_density2s4.append(density_ave2s4)


#averaging the four panels together into one because I messed up#
lst_density1 = []
lst_density2 = []
lst_density3 = []
lst_density4 = []

lst_errors1 = []
lst_errors2 = []
lst_errors3 = []
lst_errors4 = []

for i in range(13):
    lst_density1.append((lst_density1s1[i] + lst_density2s1[i] + lst_density3s1[i] + lst_density4s1[i])/4.0)
    lst_density2.append((lst_density1s2[i] + lst_density2s2[i] + lst_density3s2[i] + lst_density4s2[i])/4.0)
    lst_density3.append((lst_density1s3[i] + lst_density2s3[i] + lst_density3s3[i] + lst_density4s3[i])/4.0)
    lst_density4.append((lst_density1s4[i] + lst_density2s4[i] + lst_density3s4[i] + lst_density4s4[i])/4.0)

    lst_errors1.append((lst_errors1s1[i] + lst_errors2s1[i] + lst_errors3s1[i] + lst_errors4s1[i])/4.0)
    lst_errors2.append((lst_errors1s2[i] + lst_errors2s2[i] + lst_errors3s2[i] + lst_errors4s2[i])/4.0)
    lst_errors3.append((lst_errors1s3[i] + lst_errors2s3[i] + lst_errors3s3[i] + lst_errors4s3[i])/4.0)
    lst_errors4.append((lst_errors1s4[i] + lst_errors2s4[i] + lst_errors3s4[i] + lst_errors4s4[i])/4.0)

    
#plotting radius vs density#

pylab.errorbar(lst_r, lst_density1, yerr=lst_errors1, color='r', marker='o', markeredgecolor='none', linestyle='-', label = '9.415 < Sat. Lmass < 9.8')
pylab.errorbar(lst_r, lst_density2, yerr=lst_errors2, color='b', marker='o', markeredgecolor='none', linestyle='-', label = '9.8 < Sat. LMass < 10.3')
pylab.errorbar(lst_r, lst_density3, yerr=lst_errors3, color='g', marker='o', markeredgecolor='none', linestyle='-', label = '10.3 < Sat. LMass < 10.8')
pylab.errorbar(lst_r, lst_density4, yerr=lst_errors4, color='y', marker='o', markeredgecolor='none', linestyle='-', label = '10.8 < Sat. LMass')

pylab.suptitle('Galaxy Number Density per Aperture Radius in Four Mass Bins', fontsize=17)
pylab.xlabel("Aperture Radius (kpc)", fontsize=18)
pylab.ylabel("Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)", fontsize=18)

pylab.legend(loc=1)

pylab.xlim([10,5000])
pylab.ylim([0.00000002,0.0002])

pylab.yscale('log')
pylab.xscale('log')

pylab.ion()
pylab.show()

