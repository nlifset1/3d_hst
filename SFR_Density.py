print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math


#bring in the data#
data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_z = ascii.read("aegis_3dhst.v4.0.sfr")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0) & (data["star_flag"] == 0.0))
data_fast_flag = data_fast[idx]
data_flag = data[idx]
data_z_flag = data_z[idx]

idx2, = np.where(data_fast_flag["lsfr"] != -99)
data_fast_flagged = data_fast_flag[idx2]
data_flagged = data_flag[idx2]
data_z_flagged = data_z_flag[idx2]

errors = 0

#create a function for finding local galaxy density (projected surface density)#
def nth_nearest(gal_id, N):
    global errors
    #creating a redshift range about the chosen galaxy#
    z_un = data_z_flagged[(data_z_flagged['id'] == gal_id)]
    z_und = z_un['z']
    z = z_und[0]
    z_bin = data_fast_flagged[((data_z_flagged['z'] >= (z - 0.03)) & (data_z_flagged['z'] <= (z + 0.03)))]

    #create a list of ids of galaxies in z range#
    lst_id =[]
    for gal in z_bin:
        if gal['id'] != gal_id:
            lst_id.append(gal['id'])
    #create a list of distances from galaxy gal_id to each galaxy in z range#
    lst_dist = []
    p1 = data_flagged[(data_flagged['id'] == gal_id)]
    for gal in lst_id:
        p2 = data_flagged[(data_flagged['id'] == gal)]
        p_2 = data_z_flagged[(data_flagged['id'] == gal)]
        #convert from degrees to radians#
        ra1 = (p1['ra'])*0.0174532925
        dec1 = (p1['dec'])*0.0174532925
        dist1 = (z_un['z'])*4285
        ra2 = (p2['ra'])*0.0174532925
        dec2 = (p2['dec'])*0.0174532925
        dist2 = (p_2['z'])*4285
        #trig/math for fiding the distance between#
        x1 = (math.cos(ra1))*(math.cos(dec1))*(dist1)
        y1 = (math.sin(ra1))*(math.cos(dec1))*(dist1)
        z1 = (math.sin(ra1))*(dist1)
        x2 = (math.cos(ra2))*(math.cos(dec2))*(dist2)
        y2 = (math.sin(ra2))*(math.cos(dec2))*(dist2)
        z2 = (math.sin(ra2))*(dist2)
        #D is the distance in mpc#
        D = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        lst_dist.append(D)
    #calculate the projected surface density of galaxies around gal_id using Nth nearest#
    lst_dist.sort()
    if len(lst_dist) < (N-1):
        errors += 1
        return 0
    else:
        r_n = lst_dist[(N-1)]
        sig = N/(math.pi*(r_n**2))
        return sig

#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.4<z<2#
lst_gal_1 = []
for gal in data_fast_flagged:
    if ((gal['lmass'] >= 11.0) & ((gal['z'] >= 0.4) & (gal['z'] <= 2.0))):
        lst_gal_1.append(gal['id'])

#getting list of galaxies from previous list that also avoid edges by 0.05 degrees#
lst_gal =[]
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal)]
    if ((gal_info['ra'] < 215.25) and (gal_info['ra'] > 214.63)):
        if ((gal_info['dec'] < 53.05) and (gal_info['dec'] > 52.73)):
            lst_gal.append(gal)


#making lists for the plot#
lst_sfr = []
lst_nth = []

for gal in lst_gal:
    for item in data_fast_flagged:
        if item['id'] == gal:
            lst_sfr.append(item['lsfr'])
            lst_nth.append(nth_nearest(gal, 5))
            break

#plotting lsfr vs nth nearest#
pylab.scatter(lst_sfr, lst_nth)

pylab.suptitle('Log SFR vs Galaxy Number Density of All Redshifts', fontsize=20)
pylab.xlabel('Log Star Formation Rate', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $mpc^{-2}$)', fontsize=15)
pylab.ylim([0.0035,0.16])
pylab.yscale('log')
pylab.xlim([-3,2.5])


pylab.ion()
pylab.show()





















        
