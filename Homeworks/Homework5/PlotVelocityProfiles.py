"""
Trevor Smith, 2-28-20
ASTR400B Homework 5, Part 9

This script will do essentially the same things as PlotMassProfiles.py, with the
circular velocity profiles instead.

1. plot velocity profiles for each particle type
2. plot a total velocity profile
3. plot a visually fitted circular velocity profile based on our Hernquist profiles
from PlotMassProfiles.py
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

# import the class we created in this homework
from MassProfile import MassProfile

"""
Homework 5 Part 9:
Plot the circular velocity profile for each component of the MW up to a radius of 30 kpc.
"""

# define the range of radii on the x-axis
r = np.arange(0.1,30,0.2)*u.kpc

# use Hernquist scale lengths that we found in part 8
HernquistScales = [64.0,60.0,25.0]*u.kpc

# initialize plot
fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(15,5),
                       sharex='all',sharey='all')

# loop through galaxies (one per column)
galaxies = ['MW','M31','M33']
for i in range(len(galaxies)):
    # create MassProfile Object for this galaxy
    ProfileObj = MassProfile(galaxies[i],0)
    # get velocity profiles from particle type
    halo_v = ProfileObj.CircularVelocity(1,r)
    disk_v = ProfileObj.CircularVelocity(2,r)
    bulge_v = ProfileObj.CircularVelocity(3,r)
    # get total circular velocity profile
    total_v = ProfileObj.TotalCircularVelocity(r)

    # plot curves for velocity profiles due to each particle type
    ax[i].plot(r,halo_v,color='k',linestyle='solid',label='Dark Matter')
    ax[i].plot(r,disk_v,color='m',linestyle='solid',label='Disk Stars')
    ax[i].plot(r,bulge_v,color='r',linestyle='solid',label='Bulge Stars')
    # plot total profile
    ax[i].plot(r,total_v,color='g',linestyle='solid',label='All Matter')

    # plot circular velocity due to theoretical Hernquist profile
    # first we need the total halo mass
    ind = np.where(ProfileObj.type == 1)
    Mhalo = np.sum(ProfileObj.m[ind])*1e10*u.Msun
    HProfile = ProfileObj.HernquistVCirc(r,HernquistScales[i],Mhalo)
    ax[i].plot(r,HProfile,color='c',linestyle='dashed',
                label='Hernquist Velocities: a='+str(HernquistScales[i]))

    # set title, axis labels, and scale
    ax[i].set(title=galaxies[i]+' Circular Velocity Profile',
              xlabel='Radius (kpc)',ylabel='Velocity (km/s)',
              yscale='log')

    # show legend
    ax[i].legend()

# show plots
plt.show()
