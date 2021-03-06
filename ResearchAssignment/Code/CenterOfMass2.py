"""
Trevor Smith, 2-10-19
ASTR 400B Spring 2020
Homework 4

We would like to analyze the motion of the Local Group as a result of
gravitational forces. To do this, we must first determine the centers
of mass for the Milky Way, M31, and M33. This script will determine
the center of mass of a given galaxy at any point during the simulation.

This class assumes that the data file to be read is in the same working
directory.

"""
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy
# and simulation snapshot


    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:

        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)

        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # give astropy units to all values
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        M_total = np.sum(m)
        # xcomponent Center of mass
        Acom = np.sum(m*a)/M_total
        # ycomponent Center of mass
        Bcom = np.sum(m*b)/M_total
        # zcomponent Center of mass
        Ccom = np.sum(m*c)/M_total

        return Acom, Bcom, Ccom


    def COM_P(self,delta,VolDec):
    # Function to specifically return the center of mass position
    # iterative process to correct for mass loss and tidal disruption
    # later in the simulation
    # input:
    #   particle type (1,2,3)
    #   delta = largest permissible change in COM as we decrease the maximum
    #   radius. iterative process will stop here
    #   VolDec = fraction of previous volume we will consider each iteration
    # returns:
    #   One vector, with rows indicating:
    #   3D coordinates of the center of mass position (kpc)

        # Try a first guess at the COM position by calling COMdefine
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)

        # compute the magnitude of the COM position vector.
        RCOM = (XCOM**2 + YCOM**2 + ZCOM**2)**0.5

        # iterative process to determine the center of mass
        # change reference frame to COM frame
        # compute the difference between particle coordinates
        # and the first guess at COM position
        xNew = self.x-XCOM
        yNew = self.y-YCOM
        zNew = self.z-ZCOM
        RNEW = (xNew**2+yNew**2+zNew**2)**0.5
        # find the max 3D distance of all particles from the guessed COM
        RMAX = np.max(RNEW)

        # will re-start a reduced radius
        RMAX = RMAX/VolDec

        # pick an initial value for the change in COM position
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position
        # delta is the tolerance for the difference in the old COM and the new one.
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:
            # compute the center of mass position using
            # the particles in the reduced radius
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)

            # compute the new 3D COM position
            RCOM2 = (XCOM2**2+YCOM2**2+ZCOM2**2)**0.5

            # determine the difference between the previous center of mass position
            # and the new one.
            CHANGE = np.abs(RCOM - RCOM2)

            # uncomment the following line if you wnat to check this
            # print ("CHANGE = ", CHANGE)

            # Before loop continues, reset : RMAX, particle separations and COM
            # reduce the volume by a factor of 2 again
            RMAX = RMAX/VolDec

            # check this.
            #print ("RMAX = ", RMAX)

            # Change the frame of reference to the newly computed COM.
            # subtract the new COM
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = (xNew**2+yNew**2+zNew**2)**0.5

            # set the center of mass positions to the refined values
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        return np.around(COMP*u.kpc,2)


    def COM_V(self, COMX, COMY, COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities

        # the max distance from the center that we will use to determine the center of mass velocity
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # arguments are astropy quantities with units, so we need to give units to x,y, and z
        xV = (self.x*u.kpc)-COMX
        yV = (self.y*u.kpc)-COMY
        zV = (self.z*u.kpc)-COMZ
        RV = (xV**2+yV**2+zV**2)**0.5

        # determine the index for those particles within the max radius
        indexV = np.where(RV <= RVMAX)

        # determine the velocity and mass of those particles within the mas radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]

        # compute the center of mass velocity using those particles
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # create a vector to store the COM velocity
        # set the correct units using astropy
        # round all values
        COMV = np.array([VXCOM,VYCOM,VZCOM])
        COMV = np.around(COMV*u.km/u.s,2)

        # return the COM vector
        return COMV
