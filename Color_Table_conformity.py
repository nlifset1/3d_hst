print "start"
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import astropy.constants
from scipy import spatial
import math
from astropy.table import Table
from astropy.coordinates.sky_coordinate import SkyCoord
import random
import os
from astropy import units as u

from astropy.cosmology import WMAP9 as cosmo
os.chdir('C:\\3d_hst')

#bring in the data#
print 'getting data'
data = ascii.read('values_color.dat')
data1 = ascii.read('C:\\3d_hst/color_values.dat')


#creating a function for finding number of galaxies within a radius R (kpc)#
def Counts_q(gal_id, gal_field, z, R = 10**np.linspace(1.2,3.6,13), delta_z = 0.1, min_mass = 9.415):
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data1[data1['field'] == gal_field]

    mask = ((np.abs(data_tmp['z'] - z) <= delta_z) & (data_tmp['id'] != gal_id) & (data_tmp['lmass'] >= min_mass))
    lst_gal = data_tmp[mask]
    lst_galr = lst_gal[(((lst_gal['vj'] < 0.92) & (lst_gal['uv'] > 1.3)) | ((lst_gal['vj'] > 0.8) & (lst_gal['vj'] < 1.6) & (lst_gal['uv'] > (0.88*lst_gal['vj'] +0.49))))]
    
    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #retrieving RA and DEC data of given galaxy#
    p1 = data_tmp[(data_tmp['id'] == gal_id)]
    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(p1['ra']*u.deg, p1['dec']*u.deg)
    sc1 = SkyCoord(lst_galr['ra']*u.deg, lst_galr['dec']*u.deg)
    sc2 = SkyCoord(lst_gal['ra']*u.deg, lst_gal['dec']*u.deg)
    sep1 = sc0.separation(sc1).to(u.arcmin)
    sep2 = sc0.separation(sc2).to(u.arcmin)    
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nnr = np.empty(len(R))
    nnb = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nnr[ii] = np.sum(sep1 <= r)
        nnb[ii] = np.sum(sep2 <= r)
    lst_q = []
    for i in range(len(nnr)):
        if nnb[i] == 0:
            lst_q.append(5)
        else:
            lst_q.append(nnr[i]/nnb[i])
    return lst_q

lst_id = []
lst_field = []
lst_z = []
lst_uv = []
lst_vj = []
lst_q = []

for gal in data:
    
    lst_id.append(gal['id'])
    lst_field.append(gal['field'])
    lst_z.append(gal['z'])
    lst_uv.append(gal['uv'])
    lst_vj.append(gal['vj'])
    lst_q.append(Counts_q(gal['id'],gal['field'],gal['z']))

lst_q1 = []
lst_q2 = []
lst_q3 = []
lst_q4 = []
lst_q5 = []
lst_q6 = []
lst_q7 = []
lst_q8 = []
lst_q9 = []
lst_q10 = []
lst_q11 = []
lst_q12 = []
lst_q13 = []

for thing in lst_q:
    lst_q1.append(thing[0])
    lst_q2.append(thing[1])
    lst_q3.append(thing[2])
    lst_q4.append(thing[3])
    lst_q5.append(thing[4])
    lst_q6.append(thing[5])
    lst_q7.append(thing[6])
    lst_q8.append(thing[7])
    lst_q9.append(thing[8])
    lst_q10.append(thing[9])
    lst_q11.append(thing[10])
    lst_q12.append(thing[11])
    lst_q13.append(thing[12])

tabley = Table([lst_id,lst_field,lst_z,lst_uv,lst_vj,lst_q1,lst_q2,lst_q3,lst_q4,lst_q5,lst_q6,lst_q7,lst_q8,lst_q9,lst_q10,lst_q11,lst_q12,lst_q13], names=['id','field','z','uv','vj','q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13'])
ascii.write(tabley, 'color_tabled.dat')

