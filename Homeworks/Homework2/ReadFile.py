"""
Trevor Smith, 1-27-20
ASTR 400B Spring 2020

This function is written for HW2, so we can easily
read in our galaxy sim data
"""
import numpy as np
import astropy.units as u

def Read(filename):
    # open file in read-only mode
    file = open(filename,'r')

    # parse first line of the form: time = #
    line1 = file.readline()
    label1, value1 = line1.split()
    time = float(value1)*u.Myr

    # parse second line: num_particles = #
    line2 = file.readline()
    label2, value2 = line2.split()
    num_particles = int(value2)

    # close file
    file.close()

    # use numpy function to get the rest
    # skip first three lines
    # names = True parses the 4th line (starts with #) to get
    # column headers
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)

    # return time, # particles, and data array
    return time, num_particles, data
