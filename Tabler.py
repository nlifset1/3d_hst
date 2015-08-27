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

#create a function for finding local galaxy density (projected surface density)#
def nth_nearest(gal_id, gal_field, N):
    #creating a redshift range about the chosen galaxy#
    z_un = data_flagged[(data_flagged['id'] == gal_id) & (data_flagged['field'] == gal_field)]
    z_und = z_un['z_peak']
    z = z_und[0]
    z_bin = data_flagged[((data_flagged['z_peak'] >= (z - 0.08)) & (data_flagged['z_peak'] <= (z + 0.08)))]

    #create a list of ids of galaxies in z range and with lmass above 9.415#
    lst_id =[]
    for gal in z_bin:
        if (gal['id'] != gal_id) or (gal['field'] != gal_field):
            if gal['lmass'] >= 9.415:
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

#initializing lists for the file#
lst_id = []
lst_field = []
lst_z = []
lst_lmass = []
lst_nth = []

N = 5

#pulling necessary info for lists as well as calculating density#       
for gal in lst_gal:
    gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
    lst_id.append(gal[0])
    lst_field.append(gal[1])
    lst_z.append(gal_info['z_peak'][0])
    lst_lmass.append(gal_info['lmass'][0])
    lst_nth.append(nth_nearest(gal[0],gal[1],N))

#writing the table#
data = Table([lst_id, lst_field, lst_z, lst_lmass, lst_nth], names=['id', 'field', 'z', 'lmass', 'nth'])
ascii.write(data, 'values.dat')









