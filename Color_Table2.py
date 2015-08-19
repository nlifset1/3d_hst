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
import matplotlib.pyplot as plt
import pylab
os.chdir('C:\\3d_hst')

print 'getting data'
data = ascii.read('C:\\3d_hst/color_values.dat')

#creating a function for finding number of galaxies within a radius R (kpc)#
def Counts(gal_id, gal_field, z, R = 10**np.linspace(1.2,3.6,13), delta_z = 0.1, min_mass = 9.415):
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data[data['field'] == gal_field]

    mask = ((np.abs(data_tmp['z'] - z) <= delta_z) & (data_tmp['id'] != gal_id) & (data_tmp['lmass'] >= min_mass))
    lst_gal = data_tmp[mask]
    lst_galr = lst_gal[(((lst_gal['vj'] < 0.92) & (lst_gal['uv'] > 1.3)) | ((lst_gal['vj'] > 0.8) & (lst_gal['vj'] < 1.6) & (lst_gal['uv'] > (0.88*lst_gal['vj'] +0.49))))]
    lst_galb = lst_gal[(((lst_gal['vj'] < 0.92) & (lst_gal['uv'] < 1.3)) | ((lst_gal['vj'] > 0.8) & (lst_gal['vj'] < 1.6) & (lst_gal['uv'] < (0.88*lst_gal['vj'] +0.49))) | (lst_gal['vj']>1.5))]

    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #retrieving RA and DEC data of given galaxy#
    p1 = data_tmp[(data_tmp['id'] == gal_id)]
    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(p1['ra']*u.deg, p1['dec']*u.deg)
    sc1 = SkyCoord(lst_galr['ra']*u.deg, lst_galr['dec']*u.deg)
    sc2 = SkyCoord(lst_galb['ra']*u.deg, lst_galb['dec']*u.deg)
    sep1 = sc0.separation(sc1).to(u.arcmin)
    sep2 = sc0.separation(sc2).to(u.arcmin)    
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nnr = np.empty(len(R))
    nnb = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nnr[ii] = np.sum(sep1 <= r)
        nnb[ii] = np.sum(sep2 <= r)

    return [nnr, nnb]

#creating a function similar to Counts but for random base line#
def rand_counts(gal_field, z, R = 10**np.linspace(1.2,3.6,13), delta_z = 0.1, min_mass = 9.415):
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
    
    #switching ra and dec to degrees#
    ra1 = ra1*(180.0/math.pi)
    dec1 = dec1*(180.0/math.pi)
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data[data['field'] == gal_field]

    mask = ((np.abs(data_tmp['z'] - z) <= delta_z) & (data_tmp['lmass'] >= min_mass))
    lst_gal = data_tmp[mask]
    lst_galr = lst_gal[(((lst_gal['vj'] < 0.92) & (lst_gal['uv'] > 1.3)) | ((lst_gal['vj'] > 0.8) & (lst_gal['vj'] < 1.6) & (lst_gal['uv'] > (0.88*lst_gal['vj'] +0.49))))]
    lst_galb = lst_gal[(((lst_gal['vj'] < 0.92) & (lst_gal['uv'] < 1.3)) | ((lst_gal['vj'] > 0.8) & (lst_gal['vj'] < 1.6) & (lst_gal['uv'] < (0.88*lst_gal['vj'] +0.49))) | (lst_gal['vj']>1.5))]

    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(ra1*u.deg, dec1*u.deg)
    sc1 = SkyCoord(lst_galr['ra']*u.deg, lst_galr['dec']*u.deg)
    sc2 = SkyCoord(lst_galb['ra']*u.deg, lst_galb['dec']*u.deg)
    sep1 = sc0.separation(sc1).to(u.arcmin)
    sep2 = sc0.separation(sc2).to(u.arcmin)
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nn1 = np.empty(len(R))
    nn2 = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nn1[ii] = np.sum(sep1 <= r)
        nn2[ii] = np.sum(sep2 <= r)

    return [nn1, nn2]


def Density(gal_id, gal_field):
    gal_info = gal_info = data[(data['id'] == gal_id) & (data['field'] == gal_field)]
    lst_r = 10**np.linspace(1.2,3.6,13)
    lst_dens1 = []
    lst_dens2 = []
    lst = Counts(gal_id, gal_field, gal_info['z'])
    for i in range(len(lst_r)):
        lst_dens1.append(lst[0][i]/((lst_r[i]**2)*math.pi))
        lst_dens2.append(lst[1][i]/((lst_r[i]**2)*math.pi))
    return [lst_dens1, lst_dens2]

def Density_rand(gal_id, gal_field):
    gal_info = gal_info = data[(data['id'] == gal_id) & (data['field'] == gal_field)]
    lst_r = 10**np.linspace(1.2,3.6,13)
    lst_dens1 = []
    lst_dens2 = []
    lst1 = rand_counts(gal_field, gal_info['z'])
    lst2 = rand_counts(gal_field, gal_info['z'])
    lst = [[],[],[],[]]
    for i in range(len(lst1)):
        lst[i] = [sum(x) for x in zip(lst1[i], lst2[i])]
    for i in range(len(lst_r)):
        lst_dens1.append((lst[0][i]/((lst_r[i]**2)*math.pi))/4.0)
        lst_dens2.append((lst[1][i]/((lst_r[i]**2)*math.pi))/4.0)
    return [lst_dens1, lst_dens2]



def main():


    data = ascii.read('C:\\3d_hst/color_values.dat')
    print "Getting list of galaxies."
    #getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
    lst_gal_1 = []
    for gal in data:
        if (gal['lmass'] >= 11.0):
            if ((gal['z'] >= 0.5) and (gal['z'] <= 2.5)):
                lst_gal_1.append([gal['id'], gal['field']])

    print "There are {} galaxies in the list.".format(len(lst_gal_1))

    #EDGING#
    lst_gal = []
    for gal in lst_gal_1:
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
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
    lst_densityrr1 = []
    lst_densitybr1 = []
    lst_densityrr2 = []
    lst_densitybr2 = []
    lst_densityrr3 = []
    lst_densitybr3 = []
    lst_densityrr4 = []
    lst_densitybr4 = []
    lst_densityrr5 = []
    lst_densitybr5 = []
    lst_densityrr6 = []
    lst_densitybr6 = []
    lst_densityrr7 = []
    lst_densitybr7 = []
    lst_densityrr8 = []
    lst_densitybr8 = []
    lst_densityrr9 = []
    lst_densitybr9 = []
    lst_densityrr10 = []
    lst_densitybr10 = []
    lst_densityrr11 = []
    lst_densitybr11 = []
    lst_densityrr12 = []
    lst_densitybr12 = []
    lst_densityrr13 = []
    lst_densitybr13 = []
    
    lst_density_randrr1 = []
    lst_density_randbr1 = []
    lst_density_randrr2 = []
    lst_density_randbr2 = []
    lst_density_randrr3 = []
    lst_density_randbr3 = []
    lst_density_randrr4 = []
    lst_density_randbr4 = []
    lst_density_randrr5 = []
    lst_density_randbr5 = []
    lst_density_randrr6 = []
    lst_density_randbr6 = []
    lst_density_randrr7 = []
    lst_density_randbr7 = []
    lst_density_randrr8 = []
    lst_density_randbr8 = []
    lst_density_randrr9 = []
    lst_density_randbr9 = []
    lst_density_randrr10 = []
    lst_density_randbr10 = []
    lst_density_randrr11 = []
    lst_density_randbr11 = []
    lst_density_randrr12 = []
    lst_density_randbr12 = []
    lst_density_randrr13 = []
    lst_density_randbr13 = []
    

    lst_density = []
    lst_density_rand = []
    lst_uv =[]
    lst_vj = []


    print "Final count is {}.".format(len(lst_gal))           
    for ii, gal in enumerate(lst_gal):
        print ii
        gal_info = data[(data['id'] == gal[0]) & (data['field'] == gal[1])]
        lst_id.append(gal[0])
        lst_field.append(gal[1])
        lst_uv.append(gal_info['uv'][0])
        lst_vj.append(gal_info['vj'][0])
        lst_z.append(gal_info['z'][0])
        lst_lmass.append(gal_info['lmass'][0])
        lst_density.append(Density(gal[0],gal[1]))
        lst_density_rand.append(Density_rand(gal[0],gal[1]))

    for gal in lst_density:
        lst_densityrr1.append(gal[0][0])
        lst_densitybr1.append(gal[1][0])
        lst_densityrr2.append(gal[0][1])
        lst_densitybr2.append(gal[1][1])
        lst_densityrr3.append(gal[0][2])
        lst_densitybr3.append(gal[1][2])
        lst_densityrr4.append(gal[0][3])
        lst_densitybr4.append(gal[1][3])
        lst_densityrr5.append(gal[0][4])
        lst_densitybr5.append(gal[1][4])
        lst_densityrr6.append(gal[0][5])
        lst_densitybr6.append(gal[1][5])
        lst_densityrr7.append(gal[0][6])
        lst_densitybr7.append(gal[1][6])
        lst_densityrr8.append(gal[0][7])
        lst_densitybr8.append(gal[1][7])
        lst_densityrr9.append(gal[0][8])
        lst_densitybr9.append(gal[1][8])
        lst_densityrr10.append(gal[0][9])
        lst_densitybr10.append(gal[1][9])
        lst_densityrr11.append(gal[0][10])
        lst_densitybr11.append(gal[1][10])
        lst_densityrr12.append(gal[0][11])
        lst_densitybr12.append(gal[1][11])
        lst_densityrr13.append(gal[0][12])
        lst_densitybr13.append(gal[1][12])
    for gal in lst_density_rand:
        lst_density_randrr1.append(gal[0][0])
        lst_density_randbr1.append(gal[1][0])
        lst_density_randrr2.append(gal[0][1])
        lst_density_randbr2.append(gal[1][1])
        lst_density_randrr3.append(gal[0][2])
        lst_density_randbr3.append(gal[1][2])
        lst_density_randrr4.append(gal[0][3])
        lst_density_randbr4.append(gal[1][3])
        lst_density_randrr5.append(gal[0][4])
        lst_density_randbr5.append(gal[1][4])
        lst_density_randrr6.append(gal[0][5])
        lst_density_randbr6.append(gal[1][5])
        lst_density_randrr7.append(gal[0][6])
        lst_density_randbr7.append(gal[1][6])
        lst_density_randrr8.append(gal[0][7])
        lst_density_randbr8.append(gal[1][7])
        lst_density_randrr9.append(gal[0][8])
        lst_density_randbr9.append(gal[1][8])
        lst_density_randrr10.append(gal[0][9])
        lst_density_randbr10.append(gal[1][9])
        lst_density_randrr11.append(gal[0][10])
        lst_density_randbr11.append(gal[1][10])
        lst_density_randrr12.append(gal[0][11])
        lst_density_randbr12.append(gal[1][11])
        lst_density_randrr13.append(gal[0][12])
        lst_density_randbr13.append(gal[1][12])


    print "Gonna write table."

    data = Table([lst_id, lst_field, lst_z, lst_lmass, lst_uv, lst_vj, lst_densityrr1, lst_densitybr1, lst_densityrr2, lst_densitybr2, lst_densityrr3, lst_densitybr3, lst_densityrr4, lst_densitybr4, lst_densityrr5, lst_densitybr5, lst_densityrr6, lst_densitybr6, lst_densityrr7, lst_densitybr7, lst_densityrr8, lst_densitybr8, lst_densityrr9, lst_densitybr9, lst_densityrr10, lst_densitybr10, lst_densityrr11, lst_densitybr11, lst_densityrr12, lst_densitybr12, lst_densityrr13, lst_densitybr13, lst_density_randrr1, lst_density_randbr1, lst_density_randrr2, lst_density_randbr2, lst_density_randrr3, lst_density_randbr3, lst_density_randrr4, lst_density_randbr4, lst_density_randrr5, lst_density_randbr5, lst_density_randrr6, lst_density_randbr6, lst_density_randrr7, lst_density_randbr7, lst_density_randrr8, lst_density_randbr8, lst_density_randrr9, lst_density_randbr9, lst_density_randrr10, lst_density_randbr10, lst_density_randrr11, lst_density_randbr11, lst_density_randrr12, lst_density_randbr12, lst_density_randrr13, lst_density_randbr13], names=['id', 'field', 'z', 'lmass', 'uv', 'vj', 'densityrr1', 'densitybr1', 'densityrr2', 'densitybr2', 'densityrr3', 'densitybr3', 'densityrr4', 'densitybr4', 'densityrr5', 'densitybr5', 'densityrr6', 'densitybr6', 'densityrr7', 'densitybr7', 'densityrr8', 'densitybr8', 'densityrr9', 'densitybr9', 'densityrr10', 'densitybr10', 'densityrr11', 'densitybr11', 'densityrr12', 'densitybr12', 'densityrr13', 'densitybr13', 'randrr1', 'randbr1', 'randrr2', 'randbr2', 'randrr3', 'randbr3', 'randrr4', 'randbr4', 'randrr5', 'randbr5', 'randrr6', 'randbr6', 'randrr7', 'randbr7', 'randrr8', 'randbr8', 'randrr9', 'randbr9', 'randrr10', 'randbr10', 'randrr11', 'randbr11', 'randrr12', 'randbr12', 'randrr13', 'randbr13'])
    ascii.write(data, 'values_color.dat')
    

if __name__ == "__main__":
    main()
