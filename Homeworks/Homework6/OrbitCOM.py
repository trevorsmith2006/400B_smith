"""
Trevor Smith, 2-28-20
ASTR400B, Spring 2020

This script defines a function OrbitCOM.
This function will return the COM position and velocity vectors as a function
of time, over a given interval of snapshots.
"""

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import sys

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy = string containing galaxy name
        start = first snapshot #
        end = last snapshot #
        n = length of the intervals (in snapshots)
    returns:
        an array with columns: time, x, y, z, vx, vy, vz
        each row contains calculations for 1 snapshot
    """

    # compose the filename for output
    fileout = "Orbit_"+galaxy+".txt"

    # set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if galaxy=='M33':
        delta = 0.1
        VolDec = 4
    else:
        delta = 0.1
        VolDec = 2

    # generate the snapshot id sequence
    # check if the input is eligible
    snap_ids = np.arange(start,end,n)
    if(len(snap_ids)==0):
        sys.exit("no snapshots included in given range")

    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])

    # loop over snapshots
    for i,snap_id in enumerate(snap_ids):
        # compose the data filename
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename="VLowRes/"+"%s_"%(galaxy) + ilbl + '.txt'

        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename,2) # 2 for disk particles
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COMP = COM.COM_P(delta,VolDec)
        COMV = COM.COM_V(*COMP)

        # store the time, pos, vel in ith element of the orbit array,  without units (.value)
        # note that you can store
        # a[i] = var1, *tuple(array1)
        t = COM.time.value/1000 # units Gyr
        COMP = COMP.value
        COMV = COMV.value
        orbit[i]=t,*tuple(COMP),*tuple(COMV)

        # print snap_id to see the progress
        print(snap_id)

    # write the data to a file
    # we do this because we don't want to have to repeat this process
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first!

#OrbitCOM('M33',0,800,5)
#OrbitCOM('MW',0,800,5)
#OrbitCOM('M31',0,800,5)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

def GetPositions(galaxy):
    data = np.genfromtxt('Orbit_'+galaxy+'.txt',names=True)
    x = data['x']
    y = data['y']
    z = data['z']
    return np.array([x,y,z])
def GetVelocities(galaxy):
    data = np.genfromtxt('Orbit_'+galaxy+'.txt',names=True)
    vx = data['vx']
    vy = data['vy']
    vz = data['vz']
    return np.array([vx,vy,vz])

# function to compute the magnitude of the difference between two vectors
# You can use this function to return both the relative position and relative velocity for two
# galaxies over the entire orbit

def VectorMinus(v1,v2):
    """
    Given two arrays of 3D vectors, subtract corresponding vectors
    return an array of the magnitudes
    """
    return np.linalg.norm(v2-v1,axis=0)

# Determine the magnitude of the relative position and velocities
# of MW and M31

MWM31_dist = VectorMinus(GetPositions('MW'),GetPositions('M31'))
MWM31_vel = VectorMinus(GetVelocities('MW'),GetVelocities('M31'))

# of M33 and M31

M31M33_dist = VectorMinus(GetPositions('M31'),GetPositions('M33'))
M31M33_vel = VectorMinus(GetVelocities('M31'),GetVelocities('M33'))

# Plot the Orbit of the galaxies
#################################

# 2 rows: orbital seperation and relative orbital velocity
# 2 cols: MW vs. M31 and M31 vs. M33
fig,ax=plt.subplots(nrows=2,ncols=2,sharex='all',sharey='row')
fig.suptitle("Local Group Relative Orbits: Distance and Velocity")

# line plot: relative seperation between Milky Way and Andromeda
ax[0,0].plot(MWM31_dist)
ax[0,0].set(title="MW vs. M31",ylabel="Distance (kpc)")

# M31 and M33
ax[0,1].plot(M31M33_dist)
ax[0,1].set(title="M31 vs. M33")


# Plot the orbital velocities of the galaxies
#################################

# MW and M31
ax[1,0].plot(MWM31_vel)
ax[1,0].set(xlabel="Time (Gyr)",ylabel="Velocity (km/s)")

ax[1,1].plot(M31M33_vel)
ax[1,1].set(xlabel="Time (Gyr)")

# set log scale for y-axis (optional)
ax[0,0].set(yscale='log')
ax[1,0].set(yscale='log')

# show all plots
plt.show()
