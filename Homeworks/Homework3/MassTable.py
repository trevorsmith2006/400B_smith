"""
Trevor Smith, 2-3-20
ASTR400B, Spring 2020
Homework 3, part 3

This script will calculate the numbers discussed in part 2 of the
homework and save them in table.

table layout below

GalaxyName Mass1 Mass2 Mass3 TotalMass f_bar
MW
M31
M33

Masses 1,2, and 3 are the total masses of the three particle types
1 = Dark Matter
2 = Disk Stars
3 = Bulge Stars

TotalMass = Mass1 + Mass2 + Mass3

f_bar is the baryon fraction of a galaxy
f_bar = (Mass2+Mass3)/TotalMass
"""

import numpy as np
import astropy.units as u
from GalaxyMass import ComponentMass

# define function that returns all our column values given a data file
def CalcColumns(filepath):
    """
    This function calculates all the column values defined above.
    Inputs:
        filepath = text file that contains the galaxy data
    Outputs:
        cols = array containing all column values defined above
    """
    # mass of dark matter
    Mdark = ComponentMass(filepath,1)
    # mass of disk stars
    Mdisk = ComponentMass(filepath,2)
    # mass of bulge stars
    Mbulge = ComponentMass(filepath,3)
    # add up the total mass
    Mtotal = Mdark + Mdisk + Mbulge
    # baryonic mass fraction
    fbar = np.around((Mdisk+Mbulge)/Mtotal,3)

    #return column values in an array
    return [Mdark,Mdisk,Mbulge,Mtotal,fbar]

# save directory that contains the data files for easy reference
FilePath = '/home/tsmith/Documents/astr400b/'

# open file to write to
# we will save the data in a table using commas to seperate entries
file = open('MassTable.txt','w')

# column labels
line1='Galaxy Name, Dark Matter Mass (10^12 M_sun),'
line1+='Disk Stars Mass (10^12 M_sun),'
line1+=' Bulge Stars Mass (10^12 M_sun),'
line1+='Total Mass (10^12 M_sun),'
line1+='f_bar\n'
file.write(line1)

# Calculate Milky Way values
mwcols = CalcColumns(FilePath+'MW_000.txt')
# fill line 2 with column data returned from the function
line2 = 'MW'
for col in mwcols:
    line2+=','+str(col)
line2+='\n'

# write MW row to file
file.write(line2)

#repeat above process for M31
m31cols = CalcColumns(FilePath+'M31_000.txt')
line3 = 'M31'
for col in m31cols:
    line3 += ','+str(col)
line3+='\n'
file.write(line3)

# repeat for M33
m33cols = CalcColumns(FilePath+'M33_000.txt')
line4 = 'M33'
for col in m33cols:
    line4 += ','+str(col)
file.write(line4)

# close file
file.close()
