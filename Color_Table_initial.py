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
from operator import add
os.chdir('C:\\3d_hst')

#bring in the data#
print 'Reading in catalog...'
data_me = ascii.read("C:\\3d_hst/3dhst_big.dat")
data = ascii.read("C:\\3d_hst/3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
print 'masking'
mask = (data["use_phot"] == 1.0)
data_flagged = data_me[mask]
data_other = data[mask]
print 'there are %s galaxies' % (len(data_flagged))

mask2 = ((data_other['lmass'] >= 9.415) & (data_other['z_peak'] >= 0.4) & (data_other['z_peak'] <= 2.6))
data_color = data_flagged[mask2]
data_info = data_other[mask2]

lst_id = []
lst_field = []
lst_ra = []
lst_dec = []
lst_z = []
lst_lmass = []
lst_uv = []
lst_vj = []

for i in range(len(data_color)):
    gal = data_color[i]
    gal_info = data_info[i]
    uv = -2.5*np.log10(gal['L153']/gal['L155'])
    vj = -2.5*np.log10(gal['L155']/gal['L161'])

    lst_id.append(gal_info['id'])
    lst_field.append(gal_info['field'])
    lst_ra.append(gal_info['ra'])
    lst_dec.append(gal_info['dec'])
    lst_z.append(gal_info['z_peak'])
    lst_lmass.append(gal_info['lmass'])
    lst_uv.append(uv)
    lst_vj.append(vj)


table = Table([lst_id, lst_field, lst_ra, lst_dec, lst_z, lst_lmass, lst_uv, lst_vj], names=['id', 'field', 'ra', 'dec', 'z', 'lmass', 'uv', 'vj'])
ascii.write(table, 'color_values.dat')










