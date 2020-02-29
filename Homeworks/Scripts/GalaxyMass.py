"""
Trevor Smith, 2-3-20
ASTR 400B, Spring 2020

This script will contain a function called ComponentMass.
This function will return the total mass of any desired
galaxy component.
"""
import numpy as np
import astropy.units as u

# define our function
def ComponentMass(filename,particle_type):
    """
    Function to add up total mass of given type of matter in 
    specified galaxy data file.
    Inputs:
        filename = name of txt file with galaxy data
        particle_type = integer representing the type of matter to
        calculate mass from
    Outputs:
        total mass = total mass of the type of matter, in units of
        10^12 solar masses, rounded to 3 decimals
    """
    # read in data from file
    # numpy function will automatically organized labelled columns into
    # an array
    alldata = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)

    # save the row indices of all particles of our given type
    indices = np.where(alldata['type'] == particle_type)

    # slice an array containing the masses of these particles
    # these values are in units of 10^10 Msun
    masses = alldata['m'][indices]

    # calculate the sum of all these masses
    total_mass = np.sum(masses)

    # return this number in units of 10^12 Msun, rounded to 3 places
    # this number is already in units of 10^10 Msun
    return np.around(total_mass/1e2,3)
