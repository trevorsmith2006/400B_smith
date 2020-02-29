"""
Trevor Smith, 1-27-20
ASTR 400B Spring 2020

This funciton was written for HW2. It uses the function in ReadFile.py
to return some important information on particles in the galaxy sim.

Arguments: filename, particle type (descriptive string, see dict below)
Returns 3 data arrays for all particles of given type
1. magnitude of the distance in kpc
2. magnitude of the velocity in km/s
3. mass in solar units

"""

import numpy as np
import astropy.units as u
from ReadFile import Read

# Match description of particle types with their data value
ParticleTypes = {'Dark Matter':1.0,'Disk Stars':2.0,'Bulge Stars':3.0}

def ParticleInfo(filename,particletype):
    # read in all data
    time, n_particles, data = Read(filename)

    # save row indices of all particles of given type
    ind = np.where(data['type'] == ParticleTypes[particletype])

    # calculate and store array of distances
    # also round final magnitude to 3 decimal places
    x = data['x'][ind]
    y = data['y'][ind]
    z = data['z'][ind]
    distances = np.sqrt(x**2+y**2+z**2)*u.kpc
    distances = np.around(distances,3)

    # calculate and store velocity magnitudes
    # round to 3 places
    vx = data['vx'][ind]
    vy = data['vy'][ind]
    vz = data['vz'][ind]
    velocities = np.sqrt(vx**2+vy**2+vz**2)*(u.km/u.s)
    velocities = np.around(velocities,3)

    # store masses, in units of solar masses
    # data is stored in units of 10^10 solar masses
    masses = data['m'][ind]
    masses = masses * 10**10 * u.M_sun

    return distances, velocities, masses
