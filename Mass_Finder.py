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
        
lst_mass1 = []
lst_mass2 = []
lst_mass3 = []
lst_mass4 = []

for gal in lst_gal1:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    lst_mass1.append(gal_info['lmass'])
for gal in lst_gal2:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    lst_mass2.append(gal_info['lmass'])
for gal in lst_gal3:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    lst_mass3.append(gal_info['lmass'])
for gal in lst_gal4:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) and (data_flagged['field'] == gal[1])]
    lst_mass4.append(gal_info['lmass'])

max1 = max(lst_mass1)
max2 = max(lst_mass2)
max3 = max(lst_mass3)
max4 = max(lst_mass4)

print 'bin1: %s, bin2: %s, bin3: %s, bin4: %s' % (max1, max2, max3, max4)



                            













