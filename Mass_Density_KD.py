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

#bring in the data#
data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_z = ascii.read("aegis_3dhst.v4.0.sfr")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0) & (data_fast["lsfr"] != -99))
data_fast_flagged = data_fast[idx]
data_flagged = data[idx]
data_z_flagged = data_z[idx]


#creating a function for finding number of galaxies within a radius R (kpc)#
def Counts(gal_id, z, R):
    #making a list of galaxies in within a redshift of 0.03 of given z#
    lst_gal1 = []
    lst_gal = []
    for gal in data_z_flagged:
        if ((gal['z'] >= (z - 0.08)) and (gal['z'] <= (z + 0.08))):
            if gal['id'] != gal_id:
                lst_gal1.append(gal['id'])
    #restricting that list to galaxies above lmass of 9.415 for a 90% completeness#
    for gal in lst_gal1:
        gal_info = data_fast_flagged[(data_fast_flagged['id'] == gal)]
        if gal_info['lmass'] >= 9.415:
            lst_gal.append(gal)

    #converting radius R (kpc) to radians at given redshift z#
    kpc_per = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per = kpc_per**(-1)
    arcmin = arcmin_per*(R)
    degrees_ = arcmin/60
    degrees = degrees_.value
    radius_rad = degrees*(math.pi/180)

    #retrieving RA and DEC data (to radians) of given galaxy#
    p1 = data_flagged[(data_flagged['id'] == gal_id)]
    ra1_ = (p1['ra'])*(math.pi/180)
    ra1 = ra1_[0]
    dec1_ = (p1['dec'])*(math.pi/180)
    dec1 = dec1_[0]
    #making a list of galaxies in range of radius 'radius_rad'#
    lst_radians = []
    for gal in lst_gal:
        #pulling the necessary info of each galaxy in previous list#
        position_info = data_flagged[(data_flagged['id'] == gal)]
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

#creating a function similar to Counts but for random base line#
def rand_counts(z, R):
    #picking random location for galaxy number density#
    ra1 = random.uniform(3.746000, 3.756821)
    dec1 = random.uniform(0.920312, 0.925897)
    #making a list of galaxies in within a redshift of 0.03 of given z#
    lst_gal1 = []
    lst_gal = []
    for gal in data_z_flagged:
        if ((gal['z'] >= (z - 0.08)) and (gal['z'] <= (z + 0.08))):
            lst_gal1.append(gal['id'])
    #restricting that list to galaxies above lmass of 9.415 for a 90% completeness#
    for gal in lst_gal1:
        gal_info = data_fast_flagged[(data_fast_flagged['id'] == gal)]
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
        position_info = data_flagged[(data_flagged['id'] == gal)]
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
    if ((gal_info['z'] >= 0.5) & (gal_info['z'] <= 2.5)):
        if (gal_info['lmass'] >= 11.0):
            lst_gal_massed.append(gal)


#making lists for the plots of mass vs density, R is in MPC#
R = 20

lst_r20_density1 = []
lst_r20_density2 = []
lst_r20_density3 = []
lst_r20_density4 = []
lst_r20_sersic1 = []
lst_r20_sersic2 = []
lst_r20_sersic3 = []
lst_r20_sersic4 = []

for gal in lst_gal_massed1:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        rand_ave = 0
        for i in range(len(lst_gal_massed)):
            rand_counted += rand_counts(z, R)
        rand_ave = (float(rand_counted)/(len(lst_gal_massed)))
        lst_r20_density1.append((float(Counts(gal, z, R)) - rand_ave)/((R**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r20_lmass1.append(item['lmass'])

for gal in lst_gal_massed2:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        rand_ave = 0
        for i in range(len(lst_gal_massed)):
            rand_counted += rand_counts(z, R)
        rand_ave = (float(rand_counted)/(len(lst_gal_massed)))
        lst_r20_density2.append((float(Counts(gal, z, R)) - rand_ave)/((R**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r20_lmass2.append(item['lmass'])

for gal in lst_gal_massed3:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        rand_ave = 0
        for i in range(len(lst_gal_massed)):
            rand_counted += rand_counts(z, R)
        rand_ave = (float(rand_counted)/(len(lst_gal_massed)))
        lst_r20_density3.append((float(Counts(gal, z, R)) - rand_ave)/((R**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r20_lmass3.append(item['lmass'])

for gal in lst_gal_massed4:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        rand_ave = 0
        for i in range(len(lst_gal_massed)):
            rand_counted += rand_counts(z, R)
        rand_ave = (float(rand_counted)/(len(lst_gal_massed)))
        lst_r20_density4.append((float(Counts(gal, z, R)) - rand_ave)/((R**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r20_lmass4.append(item['lmass'])


#plotting lmass vs density#

pylab.plot(lst_r20_lmass1, lst_r20_density1, '.r-', label ='2.0 < z < 2.5')
pylab.plot(lst_r20_lmass2, lst_r20_density2, '.b-', label = '1.5 < z < 2.0')
pylab.plot(lst_r20_lmass3, lst_r20_density3, '.g-', label = '1.0 < z < 1.5')
pylab.plot(lst_r20_lmass4, lst_r20_density4, '.purple-', label = '0.5 < z < 1.0')

pylab.suptitle('Log Mass vs Log Galaxy Number Density of Four Redshift Bins', fontsize=17)
pylab.xlabel('Log Mass', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $mpc^{-2}$)', fontsize=15)
pylab.xlim([10.9, 12.0])
pylab.ylim([0.000001,0.0005])
pylab.legend(loc=1)
pylab.yscale('log')


pylab.ion()
pylab.show()

