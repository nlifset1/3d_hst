print "start"
from astropy.io import ascii
import numpy as np
import math
import random
import astropy.constants
from scipy import spatial
import os
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
os.chdir('C:\\3d_hst')

#bring in the data#
data = ascii.read("3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0))
data_flagged = data[idx]

#creating a function for finding number of galaxies within a radius R (kpc)#
def Counts(gal_id, gal_field, z, R):
    #making a list of galaxies in within a redshift of 0.08 of given z#
    lst_gal1 = []
    lst_gal = []
    for gal in data_flagged:
        if ((gal['z_peak'] >= (z - 0.08)) and (gal['z_peak'] <= (z + 0.08))):
            if (gal['id'] != gal_id) & (gal['field'] == gal_field):
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

    #retrieving RA and DEC data (to radians) of given galaxy#
    p1 = data_flagged[(data_flagged['id'] == gal_id) & (data_flagged['field'] == gal_field)]
    ra1_ = (p1['ra'])*(math.pi/180)
    ra1 = ra1_[0]
    dec1_ = (p1['dec'])*(math.pi/180)
    dec1 = dec1_[0]
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


#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z_peak'] >= 0.5) and (gal['z_peak'] <= 2.5)):
            lst_gal_1.append([gal['id'], gal['field']])

#EDGING#
lst_gal = []
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    if gal[1] == 'AEGIS':
        if (gal_info['ra'] >= 3.746/(math.pi/180)) and (gal_info['ra'] <= 3.756821/(math.pi/180)):
            if (gal_info['dec'] >= 0.920312/(math.pi/180)) and (gal_info['dec'] <= 0.925897/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'COSMOS':
        if (gal_info['ra'] >= 2.619737/(math.pi/180)) and (gal_info['ra'] <= 2.620718/(math.pi/180)):
            if (gal_info['dec'] >= 0.038741/(math.pi/180)) and (gal_info['dec'] <= 0.043811/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'GOODS-N':
        if (gal_info['ra'] >= 3.298072/(math.pi/180)) and (gal_info['ra'] <= 3.307597/(math.pi/180)):
            if (gal_info['dec'] >= 1.084787/(math.pi/180)) and (gal_info['dec'] <= 1.087936/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'GOODS-S':
        if (gal_info['ra'] >= 0.925775/(math.pi/180)) and (gal_info['ra'] <= 0.929397/(math.pi/180)):
            if (gal_info['dec'] >= -0.487098/(math.pi/180)) and (gal_info['dec'] <= -0.483591/(math.pi/180)):
                lst_gal.append(gal)
    elif gal[1] == 'UDS':
        if (gal_info['ra'] >= 0.59815/(math.pi/180)) and (gal_info['ra'] <= 0.602889/(math.pi/180)):
            if (gal_info['dec'] >= -0.091376/(math.pi/180)) and (gal_info['dec'] <= -0.090305/(math.pi/180)):
                lst_gal.append(gal)

lst_id = []
lst_field = []
lst_z = []
lst_lmass = []
lst_20 = []
lst_30 = []
lst_50 = []
lst_75 = []
lst_100 = []
lst_200 = []
lst_300 = []
lst_500 = []
lst_750 = []
lst_1000 = []

N = 5

def Density(gal_id, gal_field, R):
    gal_info = data_flagged[(data_flagged['id'] == gal_id) & (data_flagged['field'] == gal_field)]
    within = Counts(gal_id, gal_field, gal_info['z_peak'], R)
    rand_within = (1.0/3.0)*(rand_counts(gal_field, gal_info['z_peak'], R) + rand_counts(gal_field, gal_info['z_peak'], R) + rand_counts(gal_field, gal_info['z_peak'], R))
    dens = (within - rand_within)/((R**2)*math.pi)
    return dens
                
for gal in lst_gal:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_id.append(gal[0])
    lst_field.append(gal[1])
    lst_z.append(gal_info['z_peak'][0])
    lst_lmass.append(gal_info['lmass'][0])
    lst_20.append(Density(gal[0],gal[1],20))
    lst_30.append(Density(gal[0],gal[1],30))
    lst_50.append(Density(gal[0],gal[1],50))
    lst_75.append(Density(gal[0],gal[1],75))
    lst_100.append(Density(gal[0],gal[1],100))
    lst_200.append(Density(gal[0],gal[1],200))
    lst_300.append(Density(gal[0],gal[1],300))
    lst_500.append(Density(gal[0],gal[1],500))
    lst_750.append(Density(gal[0],gal[1],750))
    lst_1000.append(Density(gal[0],gal[1],1000))


data = Table([lst_id, lst_field, lst_z, lst_lmass, lst_20, lst_30, lst_50, lst_75, lst_100, lst_200, lst_300, lst_500, lst_750, lst_1000], names=['id', 'field', 'z', 'lmass', 'r20', 'r30', 'r50', 'r75', 'r100', 'r200', 'r300', 'r500', 'r750', 'r1000'])
ascii.write(data, 'values_R.dat')
