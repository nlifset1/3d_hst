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
data_s = ascii.read("aegis_3dhst.v4.1_f160w.galfit")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0) & (data_fast["lsfr"] != -99))
data_fast_flagged = data_fast[idx]
data_flagged = data[idx]
data_z_flagged = data_z[idx]
data_s_flagged = data_S[idx]


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
lst_r = [20,50,100,200]

lst_r20_density = []
lst_r50_density = []
lst_r100_density = []
lst_r200_density = []
lst_r20_sfr = []
lst_r50_sfr = []
lst_r100_sfr = []
lst_r200_sfr = []


#PLOT 1#
#R = 20 kpc#
for gal in lst_gal_massed:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        #fiding average background density#
        rand_counted = 0
        for i in range(10):
            rand_counted += rand_counts(z, 20)
        rand_ave = (float(rand_counted)/(10))
        lst_r20_density.append((float(Counts(gal, z, 20)) - rand_ave)/((20**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r20_sfr.append(item['lsfr'])
                break

#PLOT 2#
#R = 50 kpc#
for gal in lst_gal_massed:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        #fiding average background density#
        rand_counted = 0
        for i in range(10):
            rand_counted += rand_counts(z, 50)
        rand_ave = (float(rand_counted)/(10))
        lst_r50_density.append((float(Counts(gal, z, 50)) - rand_ave)/((50**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r50_sfr.append(item['lsfr'])
                break


#PLOT 3#
#R = 100 kpc#
for gal in lst_gal_massed:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        #fiding average background density#
        rand_counted = 0
        for i in range(10):
            rand_counted += rand_counts(z, 100)
        rand_ave = (float(rand_counted)/(10))
        lst_r100_density.append((float(Counts(gal, z, 100)) - rand_ave)/((100**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r100_sfr.append(item['lsfr'])
                break


#PLOT 4#
#R = 200 kpc#
for gal in lst_gal_massed:
        z_un = data_z_flagged[(data_z_flagged['id'] == gal)]
        z_und = z_un['z']
        z = z_und[0]
        #fiding average background density#
        rand_counted = 0
        for i in range(10):
            rand_counted += rand_counts(z, 200)
        rand_ave = (float(rand_counted)/(10))
        lst_r200_density.append((float(Counts(gal, z, 200)) - rand_ave)/((200**2)*math.pi))
        for item in data_fast_flagged:
            if item['id'] == gal:
                lst_r200_sfr.append(item['lsfr'])
                break


#plotting sfr vs density#

pylab.scatter(lst_sfr, lst_counts)

gs=GridSpec(2,2)
a1=pylab.subplot(gs[0,1])
a2=pylab.subplot(gs[1,1])
a3=pylab.subplot(gs[1,0])
a4=pylab.subplot(gs[1,1])

a1.plot(lst_r20_sfr, lst_r20_density, '.r-', label ='Radius 20 kpc')
a2.plot(lst_r50_sfr, lst_r50_density, '.b-', label = 'Radius 50 kpc')
a3.plot(lst_r100_sfr, lst_r100_density, '.g-', label = 'Radius 100 kpc')
a4.plot(lst_r200_sfr, lst_r200_density, '.y-', label = 'Radius 200 kpc')

pylab.suptitle('Log SFR vs Log Galaxy Number Density at Four Aperture Radii', fontsize=17)
pylab.xlabel('Log Star Formation Rate', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $mpc^{-2}$)', fontsize=15)

a1.set_xlim([-3,3])
a2.set_xlim([-3,3])
a3.set_xlim([-3,3])
a4.set_xlim([-3,3])

a1.ylim([0.0000005,0.002])
a2.ylim([0.0000005,0.002])
a3.ylim([0.0000005,0.002])
a4.ylim([0.0000005,0.002])

a1.legend(loc=1)
a2.legend(loc=1)
a3.legend(loc=1)
a4.legend(loc=1)

a1.yscale('log')
a2.yscale('log')
a3.yscale('log')
a4.yscale('log')


pylab.ion()
pylab.show()

