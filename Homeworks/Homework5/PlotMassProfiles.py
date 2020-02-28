"""
Trevor Smith, 2-28-20
ASTR400B Homework 5, Part 8

For each of our galaxies, this script will:
1. plot mass profiles for each particle type
2. plot a total mass profile
3. plot a visually fitted Hernquist mass profile to the dark matter profile
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# import the class we created in this homework
from MassProfile import MassProfile

# define the range of radii on the x-axis
r = np.arange(0.1,30,0.2)*u.kpc

# initialize plot
fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(15,5),
                       sharex='all',sharey='all')

# loop through galaxies (each column)
# also set Hernquist scale lengths for later
galaxies = ['MW','M31','M33']
HernquistScales = [64.0,60.0,25.0]*u.kpc
for i in range(len(galaxies)):
    # create MassProfile object for this galaxy
    ProfileObj = MassProfile(galaxies[i],0)
    # get component mass profiles
    halo_mass = ProfileObj.MassEnclosed(1,r)
    disk_mass = ProfileObj.MassEnclosed(2,r)
    bulge_mass = ProfileObj.MassEnclosed(3,r)
    # get total mass profile
    total_mass = ProfileObj.MassEnclosedTotal(r)

    # plot curves for component mass profiles
    ax[i].plot(r,halo_mass,color='k',linestyle='solid',label='Dark Matter')
    ax[i].plot(r,disk_mass,color='m',linestyle='solid',label='Disk Stars')
    ax[i].plot(r,bulge_mass,color='r',linestyle='solid',label='Bulge Stars')
    # plot total mass profile
    ax[i].plot(r,total_mass,color='g',linestyle='solid',label='Total Mass')

    # plot theoretical Hernquist Mass profile
    # first we need the total halo mass
    ind = np.where(ProfileObj.type == 1)
    Mhalo = np.sum(ProfileObj.m[ind])*1e10*u.Msun
    HProfile = ProfileObj.HernquistMass(r,HernquistScales[i],Mhalo)
    ax[i].plot(r,HProfile,color='c',linestyle='dashed',
               label='Hernquist Profile: a='+str(HernquistScales[i]))

    # set title, axis labels, and scale
    ax[i].set(title=galaxies[i]+' Mass Profile',
              xlabel='Radius (kpc)',ylabel='Mass Enclosed ($M_{\odot}$)',
              yscale='log')

    # show legend
    ax[i].legend()

# show final plots
plt.show()
