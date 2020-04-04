# ASTR 400B Homework 7
# 4-3-20
# Trevor Smith
# Collaborated with James Taylor

# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """

    def __init__(self,outfile):
        """
        inputs:
        outfile = string containing the filename to write orbit data to
        """
        # save output file name for later
        self.outfile = outfile

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value

        ### get the current pos/vel of M33
        M33COM = CenterOfMass('M33_000.txt',2)
        # store the position VECTOR of the M33 COM
        delta = 0.1
        VolDec = 4
        M33COMP = M33COM.COM_P(delta,VolDec)
        # store the velocity VECTOR of the M33 COM
        M33COMV = M33COM.COM_V(*M33COMP)

        # these are astropy quantities, we would like to get rid of the units
        M33COMP = M33COMP.value
        M33COMV = M33COMV.value

        ### get the current pos/vel of M31
        M31COM = CenterOfMass('M31_000.txt',2)
        # store the position VECTOR of the M31 COM
        delta = 0.1
        VolDec = 2
        M31COMP = M31COM.COM_P(delta,VolDec)
        # store the velocity VECTOR of the M31 COM
        M31COMV = M31COM.COM_V(*M31COMP)

        # these are astropy quantities, we would like to get rid of the units
        M31COMP = M31COMP.value
        M31COMV = M31COMV.value


        # store the relative position and velocity of M33 w.r.t. M31
        self.r0 = M33COMP - M31COMP
        self.v0 = M33COMV - M31COMV


        ### get the mass of each component in M31
        # disk
        self.rdisk = 5 # kpc
        self.Mdisk = ComponentMass('M31_000.txt',2)*1e12 # Msun

        # bulge
        self.rbulge = 1 # kpc
        self.Mbulge = ComponentMass('M31_000.txt',3)*1e12 # Msun

        # halo
        self.rhalo = 60 # kpc
        self.Mhalo = ComponentMass('M31_000.txt',1)*1e12 # Msun


    def HernquistAccel(self, M, r_a, r):
        """
        Calculates the acceleration due to a Hernquist mass distribution
        inputs:
        M = total mass (Msun)
        r_a = hernquist scale length (kpc)
        r = relative position vector (3D) of M33 w.r.t. M31 (kpc)
        outputs:
        3D acceleration vector (kpc/Gyr^2)
        """

        # store the magnitude of the position vector
        rmag = np.linalg.norm(r)

        # return the acceleration (3D vector) (kpc/Gyr^2)
        return -self.G*M/(rmag*(r_a+rmag)**2)*r


    def MiyamotoNagaiAccel(self,M,r_d,r):
        """
        Calculates the acceleration due to a Myamoto Nagai 1975 mass profile
        inputs:
        M = total mass (Msun)
        r_d = scale length (kpc)
        r = relative position vector (3D) (kpc)
        outputs:
        3D acceleration vector (kpc/Gyr^2)
        """
        ### Acceleration
        # calculate constants in the MN formula
        z_d = r_d/5.0
        R = (r[0]**2 + r[1]**2)**0.5
        B = r_d + (r[2]**2 + z_d**2)**0.5

        # calculate and return the acceleration vector
        # the z component is multiplied by an extra term, we will account
        # for this by multiplying by an extra array
        return -self.G*M*(R**2+B**2)**(-1.5) * r * np.array([1,1,B*(r[2]**2+z_d**2)**(-0.5)])


    def M31Accel(self,r):
        """
        Calculates the total acceleration acting on M31
        inputs:
        r = 3D position vector of M33 relative to M31
        outputs:
        3D total acceleration vector for M33 relative to M31 (kpc/Gyr^2)
        """

        ### Call the previous functions for the halo, bulge and disk
        halo_acc = self.HernquistAccel(self.Mhalo,self.rhalo,r)
        disk_acc = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)
        bulge_acc = self.HernquistAccel(self.Mbulge,self.rbulge,r)

        # return the sum
        return halo_acc+disk_acc+bulge_acc



    def LeapFrog(self, dt, r, v):
        """
        Uses the Leap Frog integration scheme to solve for the position and velocity
        of M33 at time t+dt
        inputs:
        dt = time interval (Gyr)
        r = current position vector for M33 (relative to M31) (kpc)
        v = current position vector for M33 (relative to M31) (kpc/Gyr)
        outputs:
        calculated 3D position vector of M33 at time t+dt (kpc)
        calculated 3D velocity vector of M33 at time t+dt (kpc/Gyr)
        """

        # predict the position at the next half timestep
        rhalf = r + v*dt/2

        # predict the final velocity at the next timestep using the acceleration field at the rhalf position
        vnew = v + self.M31Accel(rhalf)*dt

        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the
        # next, so we approximate it using the average expected speed over the time interval dt.
        rnew = r + 0.5*(v+vnew)*dt

        return rnew, vnew



    def OrbitIntegration(self, t0, dt, tmax):
        """
        function to caculate the projected orbit of M33 using our acceleration calculations
        numerical integration using the Leap Frog scheme
        inputs:
        t0 = initial timestep of the integration (Gyr)
        dt = time interval per integration step (Gyr)
        tmax = final timestep of the integration (Gyr)
        outputs:
        orbit = 7 column array containing 1 row per integration step
        columns: time, x, y, z, vx, vy, vz
        all units are in kpc and Gyr
        """

        # initialize the time to the input starting time
        t = t0

        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2,7))

        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]

        # initialize a counter for the orbit.
        i = 1 # since we already set the 0th values, we start the counter at 1

        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<tmax):  # as long as t has not exceeded the maximal time
            # advance the time by one timestep, dt
            t += dt
            # store the new time in the first column of the ith row
            orbit[i,0] = t

            # advance the position and velocity using the LeapFrog scheme
            r = np.array([orbit[i-1,1],orbit[i-1,2],orbit[i-1,3]])
            v = np.array([orbit[i-1,4],orbit[i-1,5],orbit[i-1,6]])
            rnew,vnew = self.LeapFrog(dt,r,v)

            # store the new position and velocity vectors into orbit[i]
            orbit[i,1:4] = rnew
            orbit[i,4:7] = vnew

            # update counter
            i+=1


        # write the data to a file
        np.savetxt(self.outfile, orbit, fmt = "%11.3f"*7, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

        # there is no return function


# Perform the integration
filename = 'AnalyticOrbit.txt'
#M33Orbit = M33AnalyticOrbit(filename)
#M33Orbit.OrbitIntegration(0,0.1,10)

# Plot Analytic Orbit
data = np.genfromtxt(filename,names=True)
fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(8,10),sharex='all')
fig.suptitle('M33 Orbit relative to M31: Analytical vs Simulated')
ax[0].plot(data['t'],np.linalg.norm([data['x'],data['y'],data['z']],axis=0),
           color='k',linestyle='--',label='analytic orbit')
#ax[0].set_xlabel('Time (Gyr)')
ax[0].set_ylabel('Seperation (kpc)')

ax[1].plot(data['t'],np.linalg.norm([data['vx'],data['vy'],data['vz']],axis=0),
           color='k',linestyle='--',label='analytic orbit')
ax[1].set_xlabel('Time(Gyr)')
ax[1].set_ylabel('Relative Velocity (kpc/Gyr)')

#Plot Simulated Orbit (results from HW 6)
M33data = np.genfromtxt('Orbit_M33.txt',names=True)
M31data = np.genfromtxt('Orbit_M31.txt',names=True)
t = M33data['t']
x = M33data['x'] - M31data['x']
y = M33data['y'] - M31data['y']
z = M33data['z'] - M31data['z']
vx = M33data['vx'] - M31data['vx']
vy = M33data['vy'] - M31data['vy']
vz = M33data['vz'] - M31data['vz']
ax[0].plot(t,np.linalg.norm([x,y,z],axis=0),
           label='simulated orbit')
ax[1].plot(t,np.linalg.norm([vx,vy,vz],axis=0),
           label='simulated orbit')
ax[0].legend()

plt.show()

"""
Homework 7 Questions

2. How do the plots compare?
The analytic orbit predicts M33 to take a very long trajectory through its
apoapsis for the second half of the simulation. The simulated orbit, however,
shows that the orbital period for M33 will dramatically increase, and it will
complete several closer orbits during this time.

3. What missing physics could make the difference?

The analytic orbit does not take into account the merger that M31 is about to
go through. M31 will dramatically increase its mass as it merges with another
galaxy of the same size. Also, the impact of the merger will probably alter the
velocity of M31 relative to M33. This could bring M33 closer in or push it farther
away depending on its position and velocity at the time.

4. The Milky Way is missing in these calculations. How might you include its effects?

That is a difficult problem to code. I would guess that we would need to alter our acceleration
calculations to account for the central potential becoming stronger as M31 gains mass. We would also
need to consider the ram pressure that M31 will experience during the merger that M33 may not experience
(or at least to the same degree). It would be difficult to model all of these things without ending up making
a simulation of our own.
"""
