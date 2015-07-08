print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math
import random

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

#creating a function for finding number of galaxies within a radius R#
def Counts(gal_id, R):
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
        dist1 = (z_un['z'])*4282.7494
        ra2 = (p2['ra'])*0.0174532925
        dec2 = (p2['dec'])*0.0174532925
        dist2 = (p_2['z'])*4282.7494
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
    #calculate the number of galaxies within the radius R#
    within = 0
    for dist in lst_dist:
        if dist <= R:
            within += 1

    return within

#creating a function similar to Counts but for random base line#
def rand_counts(R):
    #picking random location for galaxy number density#
    ra1 = random.uniform(3.746000, 3.756821)
    dec1 = random.uniform(0.920312, 0.925897)
    z = random.uniform(0.4, 2.0)
    #creating a redshift range about the chosen location#
    z_bin = data_fast_flagged[((data_z_flagged['z'] >= (z - 0.03)) & (data_z_flagged['z'] <= (z + 0.03)))]

    #create a list of ids of galaxies in z range#
    lst_id =[]
    for gal in z_bin:
        lst_id.append(gal['id'])

    #create a list of distances to each galaxy in the z bin#
    lst_dist = []
    for gal in lst_id:
        p2 = data_flagged[(data_flagged['id'] == gal)]
        p_2 = data_z_flagged[(data_flagged['id'] == gal)]
        #convert from degrees to radians#
        dist1 = (z)*4282.7494
        ra2 = (p2['ra'])*0.0174532925
        dec2 = (p2['dec'])*0.0174532925
        dist2 = (p_2['z'])*4282.7494
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

    #calculate the number of galaxies within the radius R#
    within = 0
    for dist in lst_dist:
        if dist <= R:
            within += 1

    return within
    
    return small_counts
    

#getting list of galaxies that avoid edges by 0.05 degrees#
lst_gal_edged =[]
for gal in data_flagged:
    if ((gal['ra'] < 215.25) and (gal['ra'] > 214.63)):
        if ((gal['dec'] < 53.05) and (gal['dec'] > 52.73)):
            lst_gal_edged.append(gal['id'])


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.4<z<2#
lst_gal_massed = []
for gal in lst_gal_edged:
    gal_info = data_fast_flagged[(data_z_flagged['id'] == gal)]
    if ((gal_info['z'] >= 0.4) & (gal_info['z'] <= 2.0)):
        if (gal_info['lmass'] >= 11.0):
            lst_gal_massed.append(gal)

#getting a list of all galaxies that avoid edges and aren't massive#
lst_gal_small = lst_gal_edged
for gal in lst_gal_massed:
    lst_gal_small.remove(gal)

#making lists for the plots of radius vs density, r is in MPC#
lst_r = [5,10,15,20,25,30]
lst_density = []
for r in lst_r:
    within_total = 0
    for gal in lst_gal_massed:
        within = float(Counts(gal, r) - rand_counts(r))
        other_thingy += within
    within_total = float(within_total)/float(len(lst_gal_massed))
    lst_density.append(within_total/((r**2)*math.pi))


#plotting radius vs density#

pylab.scatter(lst_r, lst_density)

pylab.suptitle('Log Mass vs Log Galaxy Number Density of All Redshifts', fontsize=23)
pylab.xlabel('Log Mass', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_gal$ $mpc^-2$)', fontsize=15)
pylab.xlim([0,30])


pylab.ion()
pylab.show()

