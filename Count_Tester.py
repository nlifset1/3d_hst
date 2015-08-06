from astropy.io import ascii
import numpy as np
import math
import random
import astropy.constants
from scipy import spatial
import os
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
os.chdir('C:\\3d_hst')

#bring in the data#
print 'Reading in catalog...'
data = ascii.read("C:\\3d_hst/3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0))
data_flagged = data[idx]


def Counts(gal_id, gal_field, z, R = 10**np.linspace(1.2,3.6,13), delta_z = 0.2, min_mass = 9.415):
    
    from astropy.coordinates.sky_coordinate import SkyCoord
    from astropy import units as u
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data_flagged[data_flagged['field'] == gal_field]

    mask = (np.abs(data_tmp['z_peak'] - z) <= delta_z) & (data_tmp['id'] != gal_id) & (data_tmp['lmass'] >= min_mass)
    lst_gal = data_tmp[mask]

    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #retrieving RA and DEC data of given galaxy#
    p1 = data_tmp[(data_tmp['id'] == gal_id)]
    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(p1['ra']*u.deg, p1['dec']*u.deg)
    sc = SkyCoord(lst_gal['ra']*u.deg, lst_gal['dec']*u.deg)
    sep = sc0.separation(sc)
    #printing random selection of seperation angles in arcminutes#
    print sep[0:80:4]
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nn = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nn[ii] = np.sum(sep <= r)
    
    return nn
