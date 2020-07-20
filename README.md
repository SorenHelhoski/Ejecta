# Ejecta
For use in analyzing the properties of ejecta from an iSALE impact:

As of now, the code is formatted to handle a spherical impactor made of one material, and up to three layers of a flat impact surface using the 2D cylindrically symmetric iSALE shock code.

Create a directory on the same level as Chicxulub, Plotting, BenchmarkData directories 

IMPORTANT: Create a text document called 'data.txt' in this new directory

Change the parameters in ejecta_search.py to the desired values

run ejecta_search.py which searches through the timesteps of 'jdata.dat' from the desired simulation

aux_library doesn't need to be run, contains Binning and KDE functions

the remaining files are different graphing programs
