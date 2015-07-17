print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math
import random
from scipy import spatial
import os
from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#making a median function easier#
def median(lst):
    return numpy.median(numpy.array(lst))

#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0))
data_flagged = data[idx]

#create a function for finding local galaxy density (projected surface density)#
def nth_nearest(gal_id, gal_field, N):
    #creating a redshift range about the chosen galaxy#
    z_un = data_flagged[(data_flagged['id'] == gal_id) and (data_flagged['field'] == gal_field)]
    z_und = z_un['z_spec']
    z = z_und[0]
    z_bin = data_fast_flagged[((data_flagged['z_spec'] >= (z - 0.08)) & (data_flagged['z_spec'] <= (z + 0.08)))]

    #create a list of ids of galaxies in z range and with lmass above 9.415#
    lst_id =[]
    for gal in z_bin:
        if (gal['id'] != gal_id) or (gal['field'] != gal_field):
            if gal['lmass'] > 9.415:
                lst_id.append([gal['id'], gal['field']])

    #finding kpc per radian ratio at given redshift z#
    kpc_arcmin = cosmo.kpc_proper_per_arcmin(z)
    kpc_degrees = kpc_arcmin*60
    kpc_radians = kpc_degrees/(math.pi/180)
    kpc_radian = kpc_radians.value
    
    #create a list of distances from galaxy gal_id to each galaxy in z range#
    lst_dist = []
    #convert from degrees to radians#
    ra1_ = (z_un['ra'])*(math.pi/180)
    ra1 = ra1_[0]
    dec1_ = (z_un['dec'])*(math.pi/180)
    dec1 = dec1_[0]
    #making a list of distances to galaxies in the z-bin in radians#
    lst_radians = []
    for gal in lst_id:
        #pulling the necessary info of each galaxy in the range#
        position_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
        ra_ = (position_info['ra'])*(math.pi/180)
        ra = ra_[0]
        dec_ = (position_info['dec'])*(math.pi/180)
        dec = dec_[0]
        #converting data to find the distance in radians to given galaxy#
        del_dec = dec - dec1
        del_ra = ra - ra1
        mean_dec = (dec + dec1)/2.0
        del_radians = math.sqrt(del_dec**2 + (del_ra*math.cos(mean_dec))**2)
        lst_radians.append(del_radians)
    lst_radians.sort()
    #finding distance to nth nearest galaxy and then calculating density from that#
    r_n_rad = lst_radians[(N-1)]
    r_n = r_n_rad*kpc_radian
    sig = N/(math.pi*(r_n**2))
    return sig

#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z_spec'] >= 0.5) & (gal['z_spec'] <= 2.5)):
            lst_gal_1.append([gal['id'], gal['field']])

#getting list of galaxies from previous list that also avoid edges by 0.05 degrees#
lst_gal =[]
upper_RA = 215.25
lower_RA = 34.27
upper_DEC = 62.33
lower_DEC = -28.0
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    if ((gal_info['ra'] < upper_RA) and (gal_info['ra'] > lower_RA)):
        if ((gal_info['dec'] < upper_DEC) and (gal_info['dec'] > lower_DEC)):
            lst_gal.append(gal)

#splitting the massed list into four redshift bins#
#bin 1: 2.0 < z < 2.5#
#bin 2: 1.5 < z < 2.0#
#bin 3: 1.0 < z < 1.5#
#bin 4: 0.5 < z < 2.0#
lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []
for gal in lst_gal:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    if gal_info['z_spec'] >= 2.0:
        lst_gal1.append(gal)
    elif gal_info['z_spec'] >= 1.5:
        lst_gal2.append(gal)
    elif gal_info['z_spec'] >= 1.0:
        lst_gal3.append(gal)
    elif gal_info['z_spec'] >= 0.5:
        lst_gal4.append(gal)

#making lists for the plot#
N = 5

#PLOT1#
#BASIC PLOT#
lst_mass1 = []
lst_nth1 = []

for gal in lst_gal1:
    for item in data_flagged:
        lst_mass1.append(item['lmass'])
        lst_nth1.append(nth_nearest(gal[0], gal[1], N))
        break

#MEDIAN POINTS#
#initializing lists for each of the 10 lmass bins#
mass_total1_1 = []
mass_total1_2 = []
mass_total1_3 = []
mass_total1_4 = []
mass_total1_5 = []
mass_total1_6 = []
mass_total1_7 = []
mass_total1_8 = []
mass_total1_9 = []
mass_total1_10 = []

density_total1_1 = []
density_total1_2 = []
density_total1_3 = []
density_total1_4 = []
density_total1_5 = []
density_total1_6 = []
density_total1_7 = []
density_total1_8 = []
density_total1_9 = []
density_total1_10 = []

#sorting the already calculated lmass and density data into the 10 lists#
for i in range(len(lst_mass1)):
    if lst_mass1[i] > 13.7:
        mass_total1_1.append(lst_mass1[i])
        density_total1_1.append(lst_nth1[i])
    elif lst_mass1[i] > 13.4:
        mass_total1_2.append(lst_mass1[i])
        density_total1_2.append(lst_nth1[i])
    elif lst_mass1[i] > 13.1:
        mass_total1_3.append(lst_mass1[i])
        density_total1_3.append(lst_nth1[i])
    elif lst_mass1[i] > 12.8:
        mass_total1_4.append(lst_mass1[i])
        density_total1_4.append(lst_nth1[i])
    elif lst_mass1[i] > 12.5:
        mass_total1_5.append(lst_mass1[i])
        density_total1_5.append(lst_nth1[i])
    elif lst_mass1[i] > 12.2:
        mass_total1_6.append(lst_mass1[i])
        density_total1_6.append(lst_nth1[i])
    elif lst_mass1[i] > 11.9:
        mass_total1_7.append(lst_mass1[i])
        density_total1_7.append(lst_nth1[i])
    elif lst_mass1[i] > 11.6:
        mass_total1_8.append(lst_mass1[i])
        density_total1_8.append(lst_nth1[i])
    elif lst_mass1[i] > 11.3:
        mass_total1_9.append(lst_mass1[i])
        density_total1_9.append(lst_nth1[i])
    elif lst_mass1[i] > 11.0:
        mass_total1_10.append(lst_mass1[i])
        density_total1_10.append(lst_nth1[i])

#calculating the median lmass and density of each bin#
density_median1_1 = median(density_total1_1)
density_median1_2 = median(density_total1_2)
density_median1_3 = median(density_total1_3)
density_median1_4 = median(density_total1_4)
density_median1_5 = median(density_total1_5)
density_median1_6 = median(density_total1_6)
density_median1_7 = median(density_total1_7)
density_median1_8 = median(density_total1_8)
density_median1_9 = median(density_total1_9)
density_median1_10 = median(density_total1_10)

mass_median1_1 = median(mass_total1_1)
mass_median1_2 = median(mass_total1_2)
mass_median1_3 = median(mass_total1_3)
mass_median1_4 = median(mass_total1_4)
mass_median1_5 = median(mass_total1_5)
mass_median1_6 = median(mass_total1_6)
mass_median1_7 = median(mass_total1_7)
mass_median1_8 = median(mass_total1_8)
mass_median1_9 = median(mass_total1_9)
mass_median1_10 = median(mass_total1_10)

#calculating standard deviation for the median points#
sigma1_1 = np.std(density_total1_1)
sigma1_2 = np.std(density_total1_2)
sigma1_3 = np.std(density_total1_3)
sigma1_4 = np.std(density_total1_4)
sigma1_5 = np.std(density_total1_5)
sigma1_6 = np.std(density_total1_6)
sigma1_7 = np.std(density_total1_7)
sigma1_8 = np.std(density_total1_8)
sigma1_9 = np.std(density_total1_9)
sigma1_10 = np.std(density_total1_10)



#PLOT2#
lst_mass2 = []
lst_nth2 = []

for gal in lst_gal2:
    for item in data_flagged:
        lst_mass2.append(item['lmass'])
        lst_nth2.append(nth_nearest(gal[0], gal[1], N))
        break

#MEDIAN POINTS#
#initializing lists for each of the 10 lmass bins#
mass_total2_1 = []
mass_total2_2 = []
mass_total2_3 = []
mass_total2_4 = []
mass_total2_5 = []
mass_total2_6 = []
mass_total2_7 = []
mass_total2_8 = []
mass_total2_9 = []
mass_total2_10 = []

density_total2_1 = []
density_total2_2 = []
density_total2_3 = []
density_total2_4 = []
density_total2_5 = []
density_total2_6 = []
density_total2_7 = []
density_total2_8 = []
density_total2_9 = []
density_total2_10 = []

#sorting the already calculated lmass and density data into the 10 lists#
for i in range(len(lst_mass2)):
    if lst_mass2[i] > 13.7:
        mass_total2_1.append(lst_mass2[i])
        density_total2_1.append(lst_nth2[i])
    elif lst_mass2[i] > 13.4:
        mass_total2_2.append(lst_mass2[i])
        density_total2_2.append(lst_nth2[i])
    elif lst_mass2[i] > 13.1:
        mass_total2_3.append(lst_mass2[i])
        density_total2_3.append(lst_nth2[i])
    elif lst_mass2[i] > 12.8:
        mass_total2_4.append(lst_mass2[i])
        density_total2_4.append(lst_nth2[i])
    elif lst_mass2[i] > 12.5:
        mass_total2_5.append(lst_mass2[i])
        density_total2_5.append(lst_nth2[i])
    elif lst_mass2[i] > 12.2:
        mass_total2_6.append(lst_mass2[i])
        density_total2_6.append(lst_nth2[i])
    elif lst_mass2[i] > 11.9:
        mass_total2_7.append(lst_mass2[i])
        density_total2_7.append(lst_nth2[i])
    elif lst_mass2[i] > 11.6:
        mass_total2_8.append(lst_mass2[i])
        density_total2_8.append(lst_nth2[i])
    elif lst_mass2[i] > 11.3:
        mass_total2_9.append(lst_mass2[i])
        density_total2_9.append(lst_nth2[i])
    elif lst_mass2[i] > 11.0:
        mass_total2_10.append(lst_mass2[i])
        density_total2_10.append(lst_nth2[i])

#calculating the median lmass and density of each bin#
density_median2_1 = median(density_total2_1)
density_median2_2 = median(density_total2_2)
density_median2_3 = median(density_total2_3)
density_median2_4 = median(density_total2_4)
density_median2_5 = median(density_total2_5)
density_median2_6 = median(density_total2_6)
density_median2_7 = median(density_total2_7)
density_median2_8 = median(density_total2_8)
density_median2_9 = median(density_total2_9)
density_median2_10 = median(density_total2_10)

mass_median2_1 = median(mass_total2_1)
mass_median2_2 = median(mass_total2_2)
mass_median2_3 = median(mass_total2_3)
mass_median2_4 = median(mass_total2_4)
mass_median2_5 = median(mass_total2_5)
mass_median2_6 = median(mass_total2_6)
mass_median2_7 = median(mass_total2_7)
mass_median2_8 = median(mass_total2_8)
mass_median2_9 = median(mass_total2_9)
mass_median2_10 = median(mass_total2_10)

#calculating standard deviation for the median points#
sigma2_1 = np.std(density_total2_1)
sigma2_2 = np.std(density_total2_2)
sigma2_3 = np.std(density_total2_3)
sigma2_4 = np.std(density_total2_4)
sigma2_5 = np.std(density_total2_5)
sigma2_6 = np.std(density_total2_6)
sigma2_7 = np.std(density_total2_7)
sigma2_8 = np.std(density_total2_8)
sigma2_9 = np.std(density_total2_9)
sigma2_10 = np.std(density_total2_10)


#PLOT3#
lst_mass3 = []
lst_nth3 = []

for gal in lst_gal3:
    for item in data_flagged:
        lst_mass3.append(item['lmass'])
        lst_nth1.append(nth_nearest(gal[0], gal[1], N))
        break

#MEDIAN POINTS#
#initializing lists for each of the 10 lmass bins#
mass_total1_1 = []
mass_total1_2 = []
mass_total1_3 = []
mass_total1_4 = []
mass_total1_5 = []
mass_total1_6 = []
mass_total1_7 = []
mass_total1_8 = []
mass_total1_9 = []
mass_total1_10 = []

density_total1_1 = []
density_total1_2 = []
density_total1_3 = []
density_total1_4 = []
density_total1_5 = []
density_total1_6 = []
density_total1_7 = []
density_total1_8 = []
density_total1_9 = []
density_total1_10 = []

#sorting the already calculated lmass and density data into the 10 lists#
for i in range(len(lst_mass3)):
    if lst_mass3[i] > 13.7:
        mass_total3_1.append(lst_mass3[i])
        density_total3_1.append(lst_nth3[i])
    elif lst_mass3[i] > 13.4:
        mass_total3_2.append(lst_mass3[i])
        density_total3_2.append(lst_nth3[i])
    elif lst_mass3[i] > 13.1:
        mass_total3_3.append(lst_mass3[i])
        density_total3_3.append(lst_nth3[i])
    elif lst_mass3[i] > 12.8:
        mass_total3_4.append(lst_mass3[i])
        density_total3_4.append(lst_nth3[i])
    elif lst_mass3[i] > 12.5:
        mass_total3_5.append(lst_mass3[i])
        density_total3_5.append(lst_nth3[i])
    elif lst_mass3[i] > 12.2:
        mass_total3_6.append(lst_mass3[i])
        density_total3_6.append(lst_nth3[i])
    elif lst_mass3[i] > 11.9:
        mass_total3_7.append(lst_mass3[i])
        density_total3_7.append(lst_nth3[i])
    elif lst_mass3[i] > 11.6:
        mass_total3_8.append(lst_mass3[i])
        density_total3_8.append(lst_nth3[i])
    elif lst_mass3[i] > 11.3:
        mass_total3_9.append(lst_mass3[i])
        density_total3_9.append(lst_nth3[i])
    elif lst_mass3[i] > 11.0:
        mass_total3_10.append(lst_mass3[i])
        density_total3_10.append(lst_nth3[i])

#calculating the median lmass and density of each bin#
density_median3_1 = median(density_total3_1)
density_median3_2 = median(density_total3_2)
density_median3_3 = median(density_total3_3)
density_median3_4 = median(density_total3_4)
density_median3_5 = median(density_total3_5)
density_median3_6 = median(density_total3_6)
density_median3_7 = median(density_total3_7)
density_median3_8 = median(density_total3_8)
density_median3_9 = median(density_total3_9)
density_median3_10 = median(density_total3_10)

mass_median3_1 = median(mass_total3_1)
mass_median3_2 = median(mass_total3_2)
mass_median3_3 = median(mass_total3_3)
mass_median3_4 = median(mass_total3_4)
mass_median3_5 = median(mass_total3_5)
mass_median3_6 = median(mass_total3_6)
mass_median3_7 = median(mass_total3_7)
mass_median3_8 = median(mass_total3_8)
mass_median3_9 = median(mass_total3_9)
mass_median3_10 = median(mass_total3_10)

#calculating standard deviation for the median points#
sigma3_1 = np.std(density_total3_1)
sigma3_2 = np.std(density_total3_2)
sigma3_3 = np.std(density_total3_3)
sigma3_4 = np.std(density_total3_4)
sigma3_5 = np.std(density_total3_5)
sigma3_6 = np.std(density_total3_6)
sigma3_7 = np.std(density_total3_7)
sigma3_8 = np.std(density_total3_8)
sigma3_9 = np.std(density_total3_9)
sigma3_10 = np.std(density_total3_10)


#PLOT4#
lst_mass4 = []
lst_nth4 = []

for gal in lst_gal4:
    for item in data_flagged:
        lst_mass4.append(item['lmass'])
        lst_nth4.append(nth_nearest(gal[0], gal[1], N))
        break

#MEDIAN POINTS#
#initializing lists for each of the 10 lmass bins#
mass_total1_1 = []
mass_total1_2 = []
mass_total1_3 = []
mass_total1_4 = []
mass_total1_5 = []
mass_total1_6 = []
mass_total1_7 = []
mass_total1_8 = []
mass_total1_9 = []
mass_total1_10 = []

density_total1_1 = []
density_total1_2 = []
density_total1_3 = []
density_total1_4 = []
density_total1_5 = []
density_total1_6 = []
density_total1_7 = []
density_total1_8 = []
density_total1_9 = []
density_total1_10 = []

#sorting the already calculated lmass and density data into the 10 lists#
for i in range(len(lst_mass4)):
    if lst_mass4[i] > 13.7:
        mass_total4_1.append(lst_mass4[i])
        density_total4_1.append(lst_nth4[i])
    elif lst_mass4[i] > 13.4:
        mass_total4_2.append(lst_mass4[i])
        density_total4_2.append(lst_nth4[i])
    elif lst_mass4[i] > 13.1:
        mass_total4_3.append(lst_mass4[i])
        density_total4_3.append(lst_nth4[i])
    elif lst_mass4[i] > 12.8:
        mass_total4_4.append(lst_mass4[i])
        density_total4_4.append(lst_nth4[i])
    elif lst_mass4[i] > 12.5:
        mass_total4_5.append(lst_mass4[i])
        density_total4_5.append(lst_nth4[i])
    elif lst_mass4[i] > 12.2:
        mass_total4_6.append(lst_mass4[i])
        density_total4_6.append(lst_nth4[i])
    elif lst_mass4[i] > 11.9:
        mass_total4_7.append(lst_mass4[i])
        density_total4_7.append(lst_nth4[i])
    elif lst_mass4[i] > 11.6:
        mass_total4_8.append(lst_mass4[i])
        density_total4_8.append(lst_nth4[i])
    elif lst_mass4[i] > 11.3:
        mass_total4_9.append(lst_mass4[i])
        density_total4_9.append(lst_nth4[i])
    elif lst_mass4[i] > 11.0:
        mass_total4_10.append(lst_mass4[i])
        density_total4_10.append(lst_nth4[i])

#calculating the median lmass and density of each bin#
density_median4_1 = median(density_total4_1)
density_median4_2 = median(density_total4_2)
density_median4_3 = median(density_total4_3)
density_median4_4 = median(density_total4_4)
density_median4_5 = median(density_total4_5)
density_median4_6 = median(density_total4_6)
density_median4_7 = median(density_total4_7)
density_median4_8 = median(density_total4_8)
density_median4_9 = median(density_total4_9)
density_median4_10 = median(density_total4_10)

mass_median4_1 = median(mass_total4_1)
mass_median4_2 = median(mass_total4_2)
mass_median4_3 = median(mass_total4_3)
mass_median4_4 = median(mass_total4_4)
mass_median4_5 = median(mass_total4_5)
mass_median4_6 = median(mass_total4_6)
mass_median4_7 = median(mass_total4_7)
mass_median4_8 = median(mass_total4_8)
mass_median4_9 = median(mass_total4_9)
mass_median4_10 = median(mass_total4_10)

#calculating standard deviation for the median points#
sigma4_1 = np.std(density_total4_1)
sigma4_2 = np.std(density_total4_2)
sigma4_3 = np.std(density_total4_3)
sigma4_4 = np.std(density_total4_4)
sigma4_5 = np.std(density_total4_5)
sigma4_6 = np.std(density_total4_6)
sigma4_7 = np.std(density_total4_7)
sigma4_8 = np.std(density_total4_8)
sigma4_9 = np.std(density_total4_9)
sigma4_10 = np.std(density_total4_10)


#PLOTTING#
#make a 2x2 plot of basic lmass vs density#

fig,pylab.axes = pylab.subplots(2, 2, sharex = True, sharey = True)

a1 = pylab.axes[0,0]
a2 = pylab.axes[0,1]
a3 = pylab.axes[1,0]
a4 = pylab.axes[1,1]

a1.plot(lst_mass1, lst_nth1, alpha=0.5, marker='o', markeredgecolor='none', linestyle='none', color='b')

a1.errorbar(mass_median1_1, density_median1_1, yerr=sigma1_1, color='r', ecolor='r')
a1.errorbar(mass_median1_2, density_median1_2, yerr=sigma1_2, color='r', ecolor='r')
a1.errorbar(mass_median1_3, density_median1_3, yerr=sigma1_3, color='r', ecolor='r')
a1.errorbar(mass_median1_4, density_median1_4, yerr=sigma1_4, color='r', ecolor='r')
a1.errorbar(mass_median1_5, density_median1_5, yerr=sigma1_5, color='r', ecolor='r')
a1.errorbar(mass_median1_6, density_median1_6, yerr=sigma1_6, color='r', ecolor='r')
a1.errorbar(mass_median1_7, density_median1_7, yerr=sigma1_7, color='r', ecolor='r')
a1.errorbar(mass_median1_8, density_median1_8, yerr=sigma1_8, color='r', ecolor='r')
a1.errorbar(mass_median1_9, density_median1_9, yerr=sigma1_9, color='r', ecolor='r')
a1.errorbar(mass_median1_10, density_median1_10, yerr=sigma1_10, color='r', ecolor='r')

a2.plot(lst_mass2, lst_nth2, alpha=0.5, marker='o', markeredgecolor='none', linestyle='none', color='b')

a1.errorbar(mass_median2_1, density_median2_1, yerr=sigma2_1, color='r', ecolor='r')
a1.errorbar(mass_median2_2, density_median2_2, yerr=sigma2_2, color='r', ecolor='r')
a1.errorbar(mass_median2_3, density_median2_3, yerr=sigma2_3, color='r', ecolor='r')
a1.errorbar(mass_median2_4, density_median2_4, yerr=sigma2_4, color='r', ecolor='r')
a1.errorbar(mass_median2_5, density_median2_5, yerr=sigma2_5, color='r', ecolor='r')
a1.errorbar(mass_median2_6, density_median2_6, yerr=sigma2_6, color='r', ecolor='r')
a1.errorbar(mass_median2_7, density_median2_7, yerr=sigma2_7, color='r', ecolor='r')
a1.errorbar(mass_median2_8, density_median2_8, yerr=sigma2_8, color='r', ecolor='r')
a1.errorbar(mass_median2_9, density_median2_9, yerr=sigma2_9, color='r', ecolor='r')
a1.errorbar(mass_median2_10, density_median2_10, yerr=sigma2_10, color='r', ecolor='r')

a3.plot(lst_mass3, lst_nth3, alpha=0.5, marker='o', markeredgecolor='none', linestyle='none', color='b')

a1.errorbar(mass_median3_1, density_median3_1, yerr=sigma3_1, color='r', ecolor='r')
a1.errorbar(mass_median3_2, density_median3_2, yerr=sigma3_2, color='r', ecolor='r')
a1.errorbar(mass_median3_3, density_median3_3, yerr=sigma3_3, color='r', ecolor='r')
a1.errorbar(mass_median3_4, density_median3_4, yerr=sigma3_4, color='r', ecolor='r')
a1.errorbar(mass_median3_5, density_median3_5, yerr=sigma3_5, color='r', ecolor='r')
a1.errorbar(mass_median3_6, density_median3_6, yerr=sigma3_6, color='r', ecolor='r')
a1.errorbar(mass_median3_7, density_median3_7, yerr=sigma3_7, color='r', ecolor='r')
a1.errorbar(mass_median3_8, density_median3_8, yerr=sigma3_8, color='r', ecolor='r')
a1.errorbar(mass_median3_9, density_median3_9, yerr=sigma3_9, color='r', ecolor='r')
a1.errorbar(mass_median3_10, density_median3_10, yerr=sigma3_10, color='r', ecolor='r')

a4.plot(lst_mass4, lst_nth4, alpha=0.5, marker='o', markeredgecolor='none', linestyle='none', color='b')

a1.errorbar(mass_median4_1, density_median4_1, yerr=sigma4_1, color='r', ecolor='r')
a1.errorbar(mass_median4_2, density_median4_2, yerr=sigma4_2, color='r', ecolor='r')
a1.errorbar(mass_median4_3, density_median4_3, yerr=sigma4_3, color='r', ecolor='r')
a1.errorbar(mass_median4_4, density_median4_4, yerr=sigma4_4, color='r', ecolor='r')
a1.errorbar(mass_median4_5, density_median4_5, yerr=sigma4_5, color='r', ecolor='r')
a1.errorbar(mass_median4_6, density_median4_6, yerr=sigma4_6, color='r', ecolor='r')
a1.errorbar(mass_median4_7, density_median4_7, yerr=sigma4_7, color='r', ecolor='r')
a1.errorbar(mass_median4_8, density_median4_8, yerr=sigma4_8, color='r', ecolor='r')
a1.errorbar(mass_median4_9, density_median4_9, yerr=sigma4_9, color='r', ecolor='r')
a1.errorbar(mass_median4_10, density_median4_10, yerr=sigma4_10, color='r', ecolor='r')

a1.set_title('0.5 < z < 1.0')
a2.set_title('1.0 < z < 1.5')
a3.set_xlabel('1.5 < z < 2.0')
a4.set_xlabel('2.0 < z < 2.5')

pylab.suptitle("Log Mass vs Galaxy Number Density in Four Redshift Bins", fontsize=17)
fig.text(0.44, 0.01, "Log Mass", fontsize=18)
fig.text(0.01, 0.67, "Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)", rotation = "vertical", fontsize=18)

pylab.xlim([10.96,14.1])
pylab.ylim([0.0000005,0.002])
pylab.yscale('log')

fig.subplots_adjust(wspace=0, hspace=0)
a1.xaxis.set_visible(False)
a2.xaxis.set_visible(False)
a2.yaxis.set_visible(False)
a4.yaxis.set_visible(False)





pylab.ion()
pylab.show()



      

























