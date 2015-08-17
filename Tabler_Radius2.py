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
data = ascii.read("C:\\3d_hst/3dhst_master.phot.v4.1.cat")

#flag out the bad stuff#
idx, = np.where((data["use_phot"] == 1.0))
data_flagged = data[idx]


#creating a function for finding number of galaxies within a radius R (kpc)#
def Counts(gal_id, gal_field, z, R = 10**np.linspace(1.2,3.6,13), delta_z = 0.1, min_mass = 9.415):
    
    from astropy.coordinates.sky_coordinate import SkyCoord
    from astropy import units as u
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data_flagged[data_flagged['field'] == gal_field]

    mask = ((np.abs(data_tmp['z_peak'] - z) <= delta_z) & (data_tmp['id'] != gal_id) & (data_tmp['lmass'] >= min_mass))
    lst_gal = data_tmp[mask]
    lst_gal1 = lst_gal[(lst_gal['lmass'] < 9.8)]
    lst_gal2 = lst_gal[((lst_gal['lmass'] < 10.3) & (lst_gal['lmass'] > 9.8))]
    lst_gal3 = lst_gal[((lst_gal['lmass'] < 10.8) & (lst_gal['lmass'] > 10.3))]
    lst_gal4 = lst_gal[((lst_gal['lmass'] < 11.8) & (lst_gal['lmass'] > 10.8))]

    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #retrieving RA and DEC data of given galaxy#
    p1 = data_tmp[(data_tmp['id'] == gal_id)]
    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(p1['ra']*u.deg, p1['dec']*u.deg)
    sc1 = SkyCoord(lst_gal1['ra']*u.deg, lst_gal1['dec']*u.deg)
    sc2 = SkyCoord(lst_gal2['ra']*u.deg, lst_gal2['dec']*u.deg)
    sc3 = SkyCoord(lst_gal3['ra']*u.deg, lst_gal3['dec']*u.deg)
    sc4 = SkyCoord(lst_gal4['ra']*u.deg, lst_gal4['dec']*u.deg)
    sep1 = sc0.separation(sc1).to(u.arcmin)
    sep2 = sc0.separation(sc2).to(u.arcmin)
    sep3 = sc0.separation(sc3).to(u.arcmin)
    sep4 = sc0.separation(sc4).to(u.arcmin)
    
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nn1 = np.empty(len(R))
    nn2 = np.empty(len(R))
    nn3 = np.empty(len(R))
    nn4 = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nn1[ii] = np.sum(sep1 <= r)
        nn2[ii] = np.sum(sep2 <= r)
        nn3[ii] = np.sum(sep3 <= r)
        nn4[ii] = np.sum(sep4 <= r)

    return [nn1, nn2, nn3, nn4]

#creating a function similar to Counts but for random base line#
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
    
    from astropy.coordinates.sky_coordinate import SkyCoord
    from astropy import units as u

    #switching ra and dec to degrees#
    ra1 = ra1*(180.0/math.pi)
    dec1 = dec1*(180.0/math.pi)
    
    #making a list of galaxies in within a redshift range of given z, in the selected field, and above the mass limit#
    lst_gal = []
    data_tmp = data_flagged[data_flagged['field'] == gal_field]

    mask = ((np.abs(data_tmp['z_peak'] - z) <= delta_z) & (data_tmp['lmass'] >= min_mass))
    lst_gal = data_tmp[mask]
    lst_gal1 = lst_gal[(lst_gal['lmass'] < 9.8)]
    lst_gal2 = lst_gal[((lst_gal['lmass'] < 10.3) & (lst_gal['lmass'] > 9.8))]
    lst_gal3 = lst_gal[((lst_gal['lmass'] < 10.8) & (lst_gal['lmass'] > 10.3))]
    lst_gal4 = lst_gal[((lst_gal['lmass'] < 11.8) & (lst_gal['lmass'] > 10.8))]

    #finding the various aperture radii in arcminutes based on given z#
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(z)
    arcmin_per_kpc = kpc_per_arcmin**(-1)
    arcmin = arcmin_per_kpc*(R*u.kpc)

    #calculating distance in special ANGLE measure to each galaxy in lst_gal#
    sc0 = SkyCoord(ra1*u.deg, dec1*u.deg)
    sc1 = SkyCoord(lst_gal1['ra']*u.deg, lst_gal1['dec']*u.deg)
    sc2 = SkyCoord(lst_gal2['ra']*u.deg, lst_gal2['dec']*u.deg)
    sc3 = SkyCoord(lst_gal3['ra']*u.deg, lst_gal3['dec']*u.deg)
    sc4 = SkyCoord(lst_gal4['ra']*u.deg, lst_gal4['dec']*u.deg)
    sep1 = sc0.separation(sc1).to(u.arcmin)
    sep2 = sc0.separation(sc2).to(u.arcmin)
    sep3 = sc0.separation(sc3).to(u.arcmin)
    sep4 = sc0.separation(sc4).to(u.arcmin)
    
    #finding number of "sep's" within the list 'arcmin' already created#
    nn1 = np.empty(len(R))
    nn2 = np.empty(len(R))
    nn3 = np.empty(len(R))
    nn4 = np.empty(len(R))
    for ii,r in enumerate(arcmin):
        nn1[ii] = np.sum(sep1 <= r)
        nn2[ii] = np.sum(sep2 <= r)
        nn3[ii] = np.sum(sep3 <= r)
        nn4[ii] = np.sum(sep4 <= r)

    return [nn1, nn2, nn3, nn4]


def Density(gal_id, gal_field):
    gal_info = gal_info = data_flagged[(data_flagged['id'] == gal_id) & (data_flagged['field'] == gal_field)]
    lst_r = 10**np.linspace(1.2,3.6,13)
    lst_dens1 = []
    lst_dens2 = []
    lst_dens3 = []
    lst_dens4 = []
    lst = Counts(gal_id, gal_field, gal_info['z_peak'])
    for i in range(len(lst_r)):
        lst_dens1.append(lst[0][i]/((lst_r[i]**2)*math.pi))
        lst_dens2.append(lst[1][i]/((lst_r[i]**2)*math.pi))
        lst_dens3.append(lst[2][i]/((lst_r[i]**2)*math.pi))
        lst_dens4.append(lst[3][i]/((lst_r[i]**2)*math.pi))
    return [lst_dens1, lst_dens2, lst_dens3, lst_dens4]

def Density_rand(gal_id, gal_field):
    gal_info = gal_info = data_flagged[(data_flagged['id'] == gal_id) & (data_flagged['field'] == gal_field)]
    lst_r = 10**np.linspace(1.2,3.6,13)
    lst_dens1 = []
    lst_dens2 = []
    lst_dens3 = []
    lst_dens4 = []
    lst1 = rand_counts(gal_field, gal_info['z_peak'])
    lst2 = rand_counts(gal_field, gal_info['z_peak'])
    lst3 = rand_counts(gal_field, gal_info['z_peak'])
    lst4 = rand_counts(gal_field, gal_info['z_peak'])
    lst = [[],[],[],[]]
    for i in range(len(lst1)):
        lst[i] = [sum(x) for x in zip(lst1[i], lst2[i], lst3[i], lst4[i])]
    for i in range(len(lst_r)):
        lst_dens1.append((lst[0][i]/((lst_r[i]**2)*math.pi))/4.0)
        lst_dens2.append((lst[1][i]/((lst_r[i]**2)*math.pi))/4.0)
        lst_dens3.append((lst[2][i]/((lst_r[i]**2)*math.pi))/4.0)
        lst_dens4.append((lst[3][i]/((lst_r[i]**2)*math.pi))/4.0)
    return [lst_dens1, lst_dens2, lst_dens3, lst_dens4]



def main():



    print "Getting list of galaxies."
    #getting a list of galaxies with lmass >= 11.0 and within the redshift range of 0.5<z<2.5#
    lst_gal_1 = []
    for gal in data_flagged:
        if (gal['lmass'] >= 11.0):
            if ((gal['z_peak'] >= 0.5) and (gal['z_peak'] <= 2.5)):
                lst_gal_1.append([gal['id'], gal['field']])

    print "There are {} galaxies in the list.".format(len(lst_gal_1))

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
    lst_density1r1 = []
    lst_density2r1 = []
    lst_density3r1 = []
    lst_density4r1 = []
    lst_density1r2 = []
    lst_density2r2 = []
    lst_density3r2 = []
    lst_density4r2 = []
    lst_density1r3 = []
    lst_density2r3 = []
    lst_density3r3 = []
    lst_density4r3 = []
    lst_density1r4 = []
    lst_density2r4 = []
    lst_density3r4 = []
    lst_density4r4 = []
    lst_density1r5 = []
    lst_density2r5 = []
    lst_density3r5 = []
    lst_density4r5 = []
    lst_density1r6 = []
    lst_density2r6 = []
    lst_density3r6 = []
    lst_density4r6 = []
    lst_density1r7 = []
    lst_density2r7 = []
    lst_density3r7 = []
    lst_density4r7 = []
    lst_density1r8 = []
    lst_density2r8 = []
    lst_density3r8 = []
    lst_density4r8 = []
    lst_density1r9 = []
    lst_density2r9 = []
    lst_density3r9 = []
    lst_density4r9 = []
    lst_density1r10 = []
    lst_density2r10 = []
    lst_density3r10 = []
    lst_density4r10 = []
    lst_density1r11 = []
    lst_density2r11 = []
    lst_density3r11 = []
    lst_density4r11 = []
    lst_density1r12 = []
    lst_density2r12 = []
    lst_density3r12 = []
    lst_density4r12 = []
    lst_density1r13 = []
    lst_density2r13 = []
    lst_density3r13 = []
    lst_density4r13 = []
    lst_density_rand1r1 = []
    lst_density_rand2r1 = []
    lst_density_rand3r1 = []
    lst_density_rand4r1 = []
    lst_density_rand1r2 = []
    lst_density_rand2r2 = []
    lst_density_rand3r2 = []
    lst_density_rand4r2 = []
    lst_density_rand1r3 = []
    lst_density_rand2r3 = []
    lst_density_rand3r3 = []
    lst_density_rand4r3 = []
    lst_density_rand1r4 = []
    lst_density_rand2r4 = []
    lst_density_rand3r4 = []
    lst_density_rand4r4 = []
    lst_density_rand1r5 = []
    lst_density_rand2r5 = []
    lst_density_rand3r5 = []
    lst_density_rand4r5 = []
    lst_density_rand1r6 = []
    lst_density_rand2r6 = []
    lst_density_rand3r6 = []
    lst_density_rand4r6 = []
    lst_density_rand1r7 = []
    lst_density_rand2r7 = []
    lst_density_rand3r7 = []
    lst_density_rand4r7 = []
    lst_density_rand1r8 = []
    lst_density_rand2r8 = []
    lst_density_rand3r8 = []
    lst_density_rand4r8 = []
    lst_density_rand1r9 = []
    lst_density_rand2r9 = []
    lst_density_rand3r9 = []
    lst_density_rand4r9 = []
    lst_density_rand1r10 = []
    lst_density_rand2r10 = []
    lst_density_rand3r10 = []
    lst_density_rand4r10 = []
    lst_density_rand1r11 = []
    lst_density_rand2r11 = []
    lst_density_rand3r11 = []
    lst_density_rand4r11 = []
    lst_density_rand1r12 = []
    lst_density_rand2r12 = []
    lst_density_rand3r12 = []
    lst_density_rand4r12 = []
    lst_density_rand1r13 = []
    lst_density_rand2r13 = []
    lst_density_rand3r13 = []
    lst_density_rand4r13 = []

    lst_density = []
    lst_density_rand = []


    print "Final count is {}.".format(len(lst_gal))           
    for ii, gal in enumerate(lst_gal):
        print ii
        gal_info = data_flagged[(data_flagged['id'] == gal[0]) & (data_flagged['field'] == gal[1])]
        lst_id.append(gal[0])
        lst_field.append(gal[1])
        lst_z.append(gal_info['z_peak'][0])
        lst_lmass.append(gal_info['lmass'][0])
        lst_density.append(Density(gal[0],gal[1]))
        lst_density_rand.append(Density_rand(gal[0],gal[1]))

    for gal in lst_density:
        lst_density1r1.append(gal[0][0])
        lst_density2r1.append(gal[1][0])
        lst_density3r1.append(gal[2][0])
        lst_density4r1.append(gal[3][0])
        lst_density1r2.append(gal[0][1])
        lst_density2r2.append(gal[1][1])
        lst_density3r2.append(gal[2][1])
        lst_density4r2.append(gal[3][1])
        lst_density1r3.append(gal[0][2])
        lst_density2r3.append(gal[1][2])
        lst_density3r3.append(gal[2][2])
        lst_density4r3.append(gal[3][2])
        lst_density1r4.append(gal[0][3])
        lst_density2r4.append(gal[1][3])
        lst_density3r4.append(gal[2][3])
        lst_density4r4.append(gal[3][3])
        lst_density1r5.append(gal[0][4])
        lst_density2r5.append(gal[1][4])
        lst_density3r5.append(gal[2][4])
        lst_density4r5.append(gal[3][4])
        lst_density1r6.append(gal[0][5])
        lst_density2r6.append(gal[1][5])
        lst_density3r6.append(gal[2][5])
        lst_density4r6.append(gal[3][5])
        lst_density1r7.append(gal[0][6])
        lst_density2r7.append(gal[1][6])
        lst_density3r7.append(gal[2][6])
        lst_density4r7.append(gal[3][6])
        lst_density1r8.append(gal[0][7])
        lst_density2r8.append(gal[1][7])
        lst_density3r8.append(gal[2][7])
        lst_density4r8.append(gal[3][7])
        lst_density1r9.append(gal[0][8])
        lst_density2r9.append(gal[1][8])
        lst_density3r9.append(gal[2][8])
        lst_density4r9.append(gal[3][8])
        lst_density1r10.append(gal[0][9])
        lst_density2r10.append(gal[1][9])
        lst_density3r10.append(gal[2][9])
        lst_density4r10.append(gal[3][9])
        lst_density1r11.append(gal[0][10])
        lst_density2r11.append(gal[1][10])
        lst_density3r11.append(gal[2][10])
        lst_density4r11.append(gal[3][10])
        lst_density1r12.append(gal[0][11])
        lst_density2r12.append(gal[1][11])
        lst_density3r12.append(gal[2][11])
        lst_density4r12.append(gal[3][11])
        lst_density1r13.append(gal[0][12])
        lst_density2r13.append(gal[1][12])
        lst_density3r13.append(gal[2][12])
        lst_density4r13.append(gal[3][12])
    for gal in lst_density_rand:
        lst_density_rand1r1.append(gal[0][0])
        lst_density_rand2r1.append(gal[1][0])
        lst_density_rand3r1.append(gal[2][0])
        lst_density_rand4r1.append(gal[3][0])
        lst_density_rand1r2.append(gal[0][1])
        lst_density_rand2r2.append(gal[1][1])
        lst_density_rand3r2.append(gal[2][1])
        lst_density_rand4r2.append(gal[3][1])
        lst_density_rand1r3.append(gal[0][2])
        lst_density_rand2r3.append(gal[1][2])
        lst_density_rand3r3.append(gal[2][2])
        lst_density_rand4r3.append(gal[3][2])
        lst_density_rand1r4.append(gal[0][3])
        lst_density_rand2r4.append(gal[1][3])
        lst_density_rand3r4.append(gal[2][3])
        lst_density_rand4r4.append(gal[3][3])
        lst_density_rand1r5.append(gal[0][4])
        lst_density_rand2r5.append(gal[1][4])
        lst_density_rand3r5.append(gal[2][4])
        lst_density_rand4r5.append(gal[3][4])
        lst_density_rand1r6.append(gal[0][5])
        lst_density_rand2r6.append(gal[1][5])
        lst_density_rand3r6.append(gal[2][5])
        lst_density_rand4r6.append(gal[3][5])
        lst_density_rand1r7.append(gal[0][6])
        lst_density_rand2r7.append(gal[1][6])
        lst_density_rand3r7.append(gal[2][6])
        lst_density_rand4r7.append(gal[3][6])
        lst_density_rand1r8.append(gal[0][7])
        lst_density_rand2r8.append(gal[1][7])
        lst_density_rand3r8.append(gal[2][7])
        lst_density_rand4r8.append(gal[3][7])
        lst_density_rand1r9.append(gal[0][8])
        lst_density_rand2r9.append(gal[1][8])
        lst_density_rand3r9.append(gal[2][8])
        lst_density_rand4r9.append(gal[3][8])
        lst_density_rand1r10.append(gal[0][9])
        lst_density_rand2r10.append(gal[1][9])
        lst_density_rand3r10.append(gal[2][9])
        lst_density_rand4r10.append(gal[3][9])
        lst_density_rand1r11.append(gal[0][10])
        lst_density_rand2r11.append(gal[1][10])
        lst_density_rand3r11.append(gal[2][10])
        lst_density_rand4r11.append(gal[3][10])
        lst_density_rand1r12.append(gal[0][11])
        lst_density_rand2r12.append(gal[1][11])
        lst_density_rand3r12.append(gal[2][11])
        lst_density_rand4r12.append(gal[3][11])
        lst_density_rand1r13.append(gal[0][12])
        lst_density_rand2r13.append(gal[1][12])
        lst_density_rand3r13.append(gal[2][12])
        lst_density_rand4r13.append(gal[3][12])


    print "Gonna write table."

    data = Table([lst_id, lst_field, lst_z, lst_lmass, lst_density1r1, lst_density1r2, lst_density1r3, lst_density1r4, lst_density1r5, lst_density1r6, lst_density1r7, lst_density1r8, lst_density1r9, lst_density1r10, lst_density1r11, lst_density1r12, lst_density1r13, lst_density2r1, lst_density2r2, lst_density2r3, lst_density2r4, lst_density2r5, lst_density2r6, lst_density2r7, lst_density2r8, lst_density2r9, lst_density2r10, lst_density2r11, lst_density2r12, lst_density2r13, lst_density3r1, lst_density3r2, lst_density3r3, lst_density3r4, lst_density3r5, lst_density3r6, lst_density3r7, lst_density3r8, lst_density3r9, lst_density3r10, lst_density3r11, lst_density3r12, lst_density3r13,lst_density4r1, lst_density4r2, lst_density4r3, lst_density4r4, lst_density4r5, lst_density4r6, lst_density4r7, lst_density4r8, lst_density4r9, lst_density4r10, lst_density4r11, lst_density4r12, lst_density4r13, lst_density_rand1r1, lst_density_rand1r2, lst_density_rand1r3, lst_density_rand1r4, lst_density_rand1r5, lst_density_rand1r6, lst_density_rand1r7, lst_density_rand1r8, lst_density_rand1r9, lst_density_rand1r10, lst_density_rand1r11, lst_density_rand1r12, lst_density_rand1r13, lst_density_rand2r1, lst_density_rand2r2, lst_density_rand2r3, lst_density_rand2r4, lst_density_rand2r5, lst_density_rand2r6, lst_density_rand2r7, lst_density_rand2r8, lst_density_rand2r9, lst_density_rand2r10, lst_density_rand2r11, lst_density_rand2r12, lst_density_rand2r13, lst_density_rand3r1, lst_density_rand3r2, lst_density_rand3r3, lst_density_rand3r4, lst_density_rand3r5, lst_density_rand3r6, lst_density_rand3r7, lst_density_rand3r8, lst_density_rand3r9, lst_density_rand3r10, lst_density_rand3r11, lst_density_rand3r12, lst_density_rand3r13, lst_density_rand4r1, lst_density_rand4r2, lst_density_rand4r3, lst_density_rand4r4, lst_density_rand4r5, lst_density_rand4r6, lst_density_rand4r7, lst_density_rand4r8, lst_density_rand4r9, lst_density_rand4r10, lst_density_rand4r11, lst_density_rand4r12, lst_density_rand4r13], names=['id', 'field', 'z', 'lmass', 'density1r1', 'density1r2', 'density1r3', 'density1r4', 'density1r5', 'density1r6', 'density1r7', 'density1r8', 'density1r9', 'density1r10', 'density1r11', 'density1r12', 'density1r13', 'density2r1', 'density2r2', 'density2r3', 'density2r4', 'density2r5', 'density2r6', 'density2r7', 'density2r8', 'density2r9', 'density2r10', 'density2r11', 'density2r12', 'density2r13', 'density3r1', 'density3r2', 'density3r3', 'density3r4', 'density3r5', 'density3r6', 'density3r7', 'density3r8', 'density3r9', 'density3r10', 'density3r11', 'density3r12', 'density3r13', 'density4r1', 'density4r2', 'density4r3', 'density4r4', 'density4r5', 'density4r6', 'density4r7', 'density4r8', 'density4r9', 'density4r10', 'density4r11', 'density4r12', 'density4r13', 'rand1r1', 'rand1r2', 'rand1r3', 'rand1r4', 'rand1r5', 'rand1r6', 'rand1r7', 'rand1r8', 'rand1r9', 'rand1r10', 'rand1r11', 'rand1r12', 'rand1r13', 'rand2r1', 'rand2r2', 'rand2r3', 'rand2r4', 'rand2r5', 'rand2r6', 'rand2r7', 'rand2r8', 'rand2r9', 'rand2r10', 'rand2r11', 'rand2r12', 'rand2r13', 'rand3r1', 'rand3r2', 'rand3r3', 'rand3r4', 'rand3r5', 'rand3r6', 'rand3r7', 'rand3r8', 'rand3r9', 'rand3r10', 'rand3r11', 'rand3r12', 'rand3r13', 'rand4r1', 'rand4r2', 'rand4r3', 'rand4r4', 'rand4r5', 'rand4r6', 'rand4r7', 'rand4r8', 'rand4r9', 'rand4r10', 'rand4r11', 'rand4r12', 'rand4r13'])
    ascii.write(data, 'values_R2.dat')


if __name__ == "__main__":
    main()
