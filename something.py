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
    #making a list of galaxies in that range#
    lst_gal = []
    for gal in data_fast_flagged:
        if ((gal['z'] >= (z - 0.03)) and (gal['z'] <= (z + 0.03))):
            if gal['id'] != gal_id:
                if gal['lmass'] > 9.415:
                    lst_gal.append(gal['id'])
    #pulling RA and DEC coordinates#
    p1 = data_flagged[(data_flagged['id'] == gal_id)]
    ra1_ = (p1['ra'])
    ra1 = ra1_[0]
    dec1_ = (p1['dec'])
    dec1 = dec1_[0]
    #making a list of difference in degree from chosen galaxy to every other in the range#
    lst_degree = []
    for gal in lst_gal:
        #pulling the necessary info#
        position_info = data_flagged[(data_flagged['id'] == gal)]
        ra_ = (position_info['ra'])
        ra = ra_[0]
        dec_ = (position_info['dec'])
        dec = dec_[0]
        #finding delta RA and delta DEC#
        del_ra = (ra-ra1)
        del_dec = (dec-dec1)
        #calculating change in degrees#
        del_degree = math.sqrt((del_ra*math.cos(dec1))**2 + del_dec**2)
        lst_degree.append(del_degree)
    #finding the kpc per arcminute at redshift z#
    kpc_per = cosmo.kpc_proper_per_arcmin(z)
    #flipping it to get arcminutes per kpc#
    min_per = kpc_per**(-1)
    #multiplying by desired radius (in kpc) to get desired aperture radius in arcminutes#
    mins = min_per*R
    #multiplying arcminutes by 0.25 to get radius in degrees#
    degrees = mins*0.25
    #calculating number of galaxies within that degree difference#
    within = 0
    for item in lst_degree:
        if item <= degrees.value:
            within += 1
    return within

#creating a function similar to Counts but for random base line#
def rand_counts(R):
    #picking random location for galaxy number density#
    ra1 = random.uniform(214.63, 215.25)
    dec1 = random.uniform(52.73, 53.05)
    z = random.uniform(0.5, 2.5)
    #making a list of galaxies in in a nearby redshift range#
    lst_gal = []
    for gal in data_fast_flagged:
        if ((gal['z'] >= (z - 0.03)) and (gal['z'] <= (z + 0.03))):
                if gal['lmass'] > 9.415:
                    lst_gal.append(gal['id'])
    #making a list of degree differences between chosen galaxy and every other in the range#
    lst_degree = []
    for gal in lst_gal:
        #pulling the necessary info#
        position_info = data_flagged[(data_flagged['id'] == gal)]
        ra_ = (position_info['ra'])
        ra = ra_[0]
        dec_ = (position_info['dec'])
        dec = dec_[0]
        #calculating delta RA and delta DEC#
        del_ra = (ra-ra1)
        del_dec = (dec-dec1)
        #calculating degree difference#
        del_degree = math.sqrt((del_ra*math.cos(dec1))**2 + del_dec**2)
        lst_degree.append(del_degree)
    #finding the kpc per arcminute at redshift z#
    kpc_per = cosmo.kpc_proper_per_arcmin(z)
    #flipping it to get arcminutes per kpc#
    min_per = kpc_per**(-1)
    #multiplying by desired radius (in kpc) to get desired aperture radius in arcminutes#
    mins = min_per*R
    #multiplying arcminutes by 0.25 to get radius in degrees#
    degrees = mins*0.25
    #finding number of galaxies within range#
    within = 0
    for item in lst_degree:
        if item <= degrees.value:
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

#getting a list of all galaxies that avoid edges and aren't massive#
lst_gal_small = lst_gal_edged
for gal in lst_gal_massed:
    lst_gal_small.remove(gal)

#making lists for the plots of mass vs density, R is in KPC#
R = 10000
rand_counted = 0
for i in range(len(lst_gal_massed)):
    rand_counted += rand_counts(R)
rand_counted = float(rand_counted/len(lst_gal_massed))
lst_mass =[]
lst_counts =[]
for gal in lst_gal_massed:
    gal_counted = Counts(gal, R)
    lst_counts.append((float(gal_counted) - rand_counted)/((R**2)*math.pi))
    for item in data_fast_flagged:
        if item['id'] == gal:
            lst_mass.append(item['lmass'])


#plotting mass vs density#

pylab.scatter(lst_mass, lst_counts)

pylab.suptitle('Log Mass vs Log Galaxy Number Density of All Redshifts', fontsize=20)
pylab.xlabel('Log Mass', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.xlim([10.96,11.6])
pylab.ylim([0.000000005,0.0000005])
pylab.yscale('log')


pylab.ion()
pylab.show()

