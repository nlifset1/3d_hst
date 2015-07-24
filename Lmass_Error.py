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
stuff = ascii.read('values_R.dat')

#flag out the bad stuff#
data_flagged = data[(data["use_phot"] == 1.0)]

#creating a function similar to Counts but for random base line#
def rand_counts(gal_field, z, R):
    #picking random location for galaxy number density#
    if gal_field == 'AEGIS':
        ra1 = random.uniform(3.746000, 3.756821)
        dec1 = random.uniform(0.920312, 0.925897)
    elif gal_field == 'COSMOS':
        ra1 = random.uniform(2.619737, 2.620718)
        dec1 = random.uniform(0.038741, 0.043811)
    elif gal_field == 'GOODS-N':
        ra1 = random.uniform(3.298072, 3.307597)
        dec1 = random.uniform(1.084787, 1.087936)
    elif gal_field == 'GOODS-S':
        ra1 = random.uniform(0.925775, 0.929397)
        dec1 = random.uniform(-0.487098, -0.483591)
    elif gal_field == 'UDS':
        ra1 = random.uniform(0.59815, 0.602889)
        dec1 = random.uniform(-0.091376, -0.090305)
    
    #making a list of galaxies in within a redshift of 0.03 of given z#
    lst_gal1 = []
    lst_gal = []
    for gal in data_flagged:
        if gal['field'] == gal_field:
            if ((gal['z_peak'] >= (z - 0.08)) and (gal['z_peak'] <= (z + 0.08))):
                lst_gal1.append([gal['id'], gal['field']])
    #restricting that list to galaxies above lmass of 9.415 for a 90% completeness#
    for gal in lst_gal1:
        gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
        if gal_info['lmass'] >= 9.415:
            lst_gal.append(gal)

    #converting radius R (kpc) to radians at given redshift z#
    kpc_per = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per = kpc_per**(-1)
    arcmin = arcmin_per*(R)
    degrees_ = arcmin/60
    degrees = degrees_.value
    radius_rad = degrees*(math.pi/180)

    #making a list of galaxies in range of radius 'radius_rad'#
    lst_radians = []
    for gal in lst_gal:
        #pulling the necessary info of each galaxy in previous list#
        position_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
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
    #finding number of distances in lst_radians that are within calculated radius_rad#
    within = 0
    for dist in lst_radians:
        if dist <= radius_rad:
            within += 1
    return within

lst_errors1 = []
lst_errors2 = []
lst_errors3 = []
lst_errors4 = []

lst_radius = [20,30,50,75,100,200,300,500,750,1000]

lst_gal1 = []
lst_gal2 = []
lst_gal3 = []
lst_gal4 = []

for gal in stuff:
    if gal['lmass'] <= 11.15:
        lst_gal1.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.3:
        lst_gal2.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.45:
        lst_gal3.append([gal['id'],gal['field']])
    elif gal['lmass'] <= 11.6:
        lst_gal4.append([gal['id'],gal['field']])

for R in lst_radius:

    lst_rand1 = []
    for gal in lst_gal1:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand.append(rand_counts(gal[1], gal_info['z'], R))
    lst_errors1.append(np.std(lst_rand1))

    lst_rand2 = []
    for gal in lst_gal2:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand.append(rand_counts(gal[1], gal_info['z'], R))
    lst_errors2.append(np.std(lst_rand2))

    lst_rand3 = []
    for gal in lst_gal3:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand.append(rand_counts(gal[1], gal_info['z'], R))
    lst_errors3.append(np.std(lst_rand3))

    lst_rand4 = []
    for gal in lst_gal4:
        gal_info = stuff[(stuff['id'] == gal[0]) & (stuff['field'] == gal[1])]
        lst_rand.append(rand_counts(gal[1], gal_info['z'], R))
    lst_errors4.append(np.std(lst_rand4))

data = Table([lst_errors1,lst_errors2,lst_errors3,lst_errors4], names=['std1','std2','std3','std4'])
ascii.write(data, 'lmass_errors.dat')



