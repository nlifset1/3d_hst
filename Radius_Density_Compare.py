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

#making lists for the plots of radius vs density, r is in KPC#
lst_r = [20,50,100,200]
lst_density = []
lst_rand = []
for r in lst_r:
    density_total = 0
    density_rand_total = 0
    for gal in lst_gal_massed:
        #finding number density of each galaxy at given radius and averaging them#
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        within = float(Counts(gal, z, r))
        within_rand = float(rand_counts(z, r)+rand_counts(z, r))/2.0
        density = within/((r**2)*math.pi)
        density_rand = within_rand/((r**2.0)*math.pi)
        density_total += density
        density_rand_total += density_rand
    density_ave = float(density_total)/len(lst_gal_massed)
    density_rand_ave = float(density_rand_total)/len(lst_gal_massed)
    lst_density.append(density_ave)
    lst_rand.append(density_rand_ave)


#plotting radius vs density#

pylab.plot(lst_r, lst_density, marker = 'o', markeredgecolor='none', color='b', linestyle='-', label = 'Selected Massive Galaxies')
pylab.plot(lst_r, lst_rand, marker = 'o', markeredgecolor='none', color='r', linestyle='-', label = 'Average Background Number Density')

pylab.suptitle('Galaxy Number Density per Aperture Radius at All Redshifts', fontsize=17)
pylab.xlabel('Aperture Radius (kpc)', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.xlim([16,210])
pylab.ylim([0.000001,0.0002])
pylab.yscale('log')
pylab.xscale('log')
pylab.legend(loc=1)


pylab.ion()
pylab.show()

