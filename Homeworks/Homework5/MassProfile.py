"""
Trevor Smith, 2-16-20
ASTR400B Homework 5

This script will define an importable class MassProfile.

Collaborated with Michael Klein and James Taylor.
"""
import numpy as np
import astropy.units as u

# import our functions from previous homeworks
from ReadFile import Read
from CenterOfMass import CenterOfMass

# import constants
import astropy.constants as const

class MassProfile:
    def __init__(self, galaxy, snap):
        """
        Inputs:
            galaxy = string containing abbreviated galaxy name
            examples are MW, M31,M33
            snap = snap number defining the timestep
        """
        # reconstruct data filename from args
        # add a string of the filenumber to the value "000"
        snapstring = '000' + str(snap)
        # remove all but the last 3 digits
        snapstring = snapstring[-3:]
        self.filename = "%s_"%(galaxy) + snapstring + ".txt"

        # read in data
        time, total, data = Read(self.filename)
        # save and give x, y, z appropriate units
        self.x = data['x']*u.kpc
        self.y = data['y']*u.kpc
        self.z = data['z']*u.kpc

        # save mass (in 1e10 Msun, but no astropy units)
        self.m = data['m']

        # save particle types
        self.type = data['type']

        # save galaxy name
        self.gname = galaxy

    def MassEnclosed(self, ptype, r):
        """
        Calculates the enclosed mass for a specific particle type.
        Inputs:
            ptype = particle type
            r = array of radii to calculate enclosed mass (u.kpc)
        Outputs:
            array of enclosed masses (u.Msun)
        """
        # M33 does not have any bulge stars
        # check for this case before any calculations
        if self.gname=='M33' and ptype==3:
            return np.zeros(len(r))

        # create CenterOfMass object
        # always use disk stars for CoM
        COM = CenterOfMass(self.filename, 2) # ptype=2 for disk stars

        # get COM position
        delta = 0.1
        COMX, COMY, COMZ = COM.COM_P(delta)

        # get particle positions in COM frame of reference
        xNew = self.x - COMX
        yNew = self.y - COMY
        zNew = self.z - COMZ

        # select only particles of ptype
        ind = np.where(self.type == ptype)
        xNew2 = xNew[ind]
        yNew2 = yNew[ind]
        zNew2 = zNew[ind]
        mNew = self.m[ind]

        # Calculate magnitudes
        R = (xNew2**2+yNew2**2+zNew2**2)**(0.5)
        # create array to store answers
        encmass = np.zeros(len(r))

        # loop through r to do the calculations
        for i in range(len(r)):
            # find indices of particles within this radius
            index = np.where(R < r[i])

            # add their masses and store in our array
            encmass[i] = np.sum(mNew[index])

        # return masses with proper units
        return encmass*1e10*u.Msun

    def MassEnclosedTotal(self,r):
        """
        Calculates the total mass enclosed within a given radius
        Inputs:
            r = array of radii to calculate total enclosed mass (u.kpc)
        Outputs:
            array of enclosed masses (u.Msun)
        """

        # get enclosed masses at r for all particle types
        DiskMass = self.MassEnclosed(2,r)
        DarkMass = self.MassEnclosed(1,r)
        BulgeMass = self.MassEnclosed(3,r)

        # return the sum of these at each r
        return DarkMass+DiskMass+BulgeMass

    def HernquistMass(self,r,a,Mhalo):
        """
        Computes the mass enclosed within a given radius
        according to the theoretical Hernquist mass profile.
        Inputs:
            r = radius to calculate enclosed mass within (u.kpc)
            a = scale length for Hernquist profile (u.kpc)
            Mhalo = total dark matter halo mass (u.Msun)
        Outputs:
            enclosed mass (u.Msun)
        """
        return Mhalo*r**2/(a+r)**2

    def CircularVelocity(self,ptype,r):
        """
        Assuming spherical symmetry, this function calculates the circular speed
        due to a specific type of particle at a given radius.
        Inputs:
            ptype = particle type
            r = array of radii (u.kpc)
        Outputs:
            array of circular speeds (u.km/u.s)
        """
        # define G in proper units for this calculation
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # plug into equation, answer should be in km/s
        # round to 2 decimals
        return np.around(np.sqrt(G*self.MassEnclosed(ptype,r)/r),2)

    def TotalCircularVelocity(self,r):
        """
        Assuming spherical symmetry, this function calculates the circular speed
        as a result of all particles at a given radius.
        Inputs:
            r = array of radii (u.kpc)
        Outputs:
            array of circular speeds (u.km/u.s)
        """
        # define G in proper units for this calculation
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # plug into equation, answer should be in km/s
        # round to 2 decimals
        return np.around(np.sqrt(G*self.MassEnclosedTotal(r)/r),2)

    def HernquistVCirc(self,r,a,Mhalo):
        """
        Assuming spherical symmetry, calculates the circular speed at a given
        radius. This function assumes a Hernquist mass profile.
        Inputs:
            r = array of radii (u.kpc)
            a = scale factor for Hernquist profile (u.kpc)
            Mhalo = mass of dark matter halo (u.Msun)
        Outputs:
            array of velocities (u.km/u.s)
        """
        # define G in useful units
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # plug into equation, round to 2 decimal places
        return np.around(np.sqrt(G*self.HernquistMass(r,a,Mhalo)/r),2)
