"""
Trevor Smith, 1-27-20
ASTR 400B, Spring 2020

This script will test the functions ParticleInfo in 
ParticleProperties.py and Read in ReadFile.py.
We will also attempt to change units using Astropy built-in
functions.
"""
import numpy as np
import astropy.units as u
from ParticleProperties import ParticleInfo

# We will use this linked file to test the code
filepath = '/home/tsmith/Documents/astr400b/MW_000.txt'

# call constructed function to retrieve data
dist, vel, mass = ParticleInfo(filepath,'Disk Stars')

# print information about 100th disk particle
# index 99
print("Printing data for 100th disk particle:")
print("3D distance = ",dist[99])
print("3D velocity = ",vel[99])
print("Mass = ",mass[99])

# use Quantitiy.to() to convert distance to LYs
print("Converting distance to LYs")
newdist = dist[99].to(u.lyr)
print("3D distance = ",newdist)
