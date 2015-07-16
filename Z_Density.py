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

#create a function for finding local galaxy density (projected surface density)#
def nth_nearest(gal_id, N):
    #creating a redshift range about the chosen galaxy#
    z_un = data_z_flagged[(data_z_flagged['id'] == gal_id)]
    z_und = z_un['z']
    z = z_und[0]
    z_bin = data_fast_flagged[((data_z_flagged['z'] >= (z - 0.08)) & (data_z_flagged['z'] <= (z + 0.08)))]

    #create a list of ids of galaxies in z range and with lmass above 9.415#
    lst_id =[]
    for gal in z_bin:
        if gal['id'] != gal_id:
            if gal['lmass'] > 9.415:
                lst_id.append(gal['id'])

    #finding kpc per radian ratio at given redshift z#
    kpc_arcmin = cosmo.kpc_proper_per_arcmin(z)
    kpc_degrees = kpc_arcmin*60
    kpc_radians = kpc_degrees/(math.pi/180)
    kpc_radian = kpc_radians.value
    
    #create a list of distances from galaxy gal_id to each galaxy in z range#
    lst_dist = []
    p1 = data_flagged[(data_flagged['id'] == gal_id)]
    #convert from degrees to radians#
    ra1_ = (p1['ra'])*(math.pi/180)
    ra1 = ra1_[0]
    dec1_ = (p1['dec'])*(math.pi/180)
    dec1 = dec1_[0]
    #making a list of distances to galaxies in the z-bin in radians#
    lst_radians = []
    for gal in lst_gal:
        #pulling the necessary info of each galaxy in the range#
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
    lst_radians.sort()
    #finding distance to nth nearest galaxy and then calculating density from that#
    r_n_rad = lst_radians[(N-1)]
    r_n = r_n_rad*kpc_radian
    sig = N/(math.pi*(r_n**2))
    return sig

#getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
lst_gal_1 = []
for gal in data_fast_flagged:
    if (gal['lmass'] >= 11.0):
        if ((gal['z'] >= 0.5) & (gal['z'] <= 2.5)):
            lst_gal_1.append(gal['id'])

#getting list of galaxies from previous list that also avoid edges by 0.05 degrees#
lst_gal =[]
for gal in lst_gal_1:
    gal_info = data_flagged[(data_flagged['id'] == gal)]
    if ((gal_info['ra'] < 215.25) and (gal_info['ra'] > 214.63)):
        if ((gal_info['dec'] < 53.05) and (gal_info['dec'] > 52.73)):
            lst_gal.append(gal)


#making lists for the plot#
lst_z = []
lst_nth = []
N = 5

for gal in lst_gal:
    for item in data_z_flagged:
        if item['id'] == gal:
            lst_z.append(item['z'])
            lst_nth.append(nth_nearest(gal, N))
            break
        

#plotting redshift vs nth nearest#
pylab.scatter(lst_z, lst_nth, alpha=0.7)

pylab.suptitle('Redshift vs Galaxy Number Density', fontsize=20)
pylab.xlabel('Redshift', fontsize=16)
pylab.ylabel('Log Galaxy Number Density ($N_{gal}$ $kpc^{-2}$)', fontsize=15)
pylab.xlim([0.4,2.6])
pylab.ylim([0.0000005,0.002])
pylab.yscale('log')


pylab.ion()
pylab.show()
