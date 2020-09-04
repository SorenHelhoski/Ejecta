# Ejecta
For use in analyzing the properties of ejecta from an iSALE impact:

As of now, the code is formatted to handle a spherical impactor made of one material, and up to three layers of a flat impact surface using the 2D cylindrically symmetric iSALE shock code. Tracers must be enabled, and set to -1 spacing in all directions.

Create a directory on the same level as (ImpactDataDir), Plotting, BenchmarkData, and eos directories 

Change the parameters in Ejecta.inp to the desired values

run ejecta_search.py which searches through the timesteps of 'jdata.dat' from the desired simulation

Rough description of method: tracers passing above y=0 are appended to an ejecta array containing their starting x and y coord, launch angle and velocity, landing location, material, and peak temperature and pressures. (Pressure is taken as the max over the entire simulation, temperature is taken as the maximum up until launch)

aux_library doesn't need to be run, contains Binning and KDE functions

Launch.py : mass density dist, velocity & angle dist 

Landing.py : landing loc dist, pres & temp dist, cumulative ejecta 

Panels.py : KDEs and histograms of the peak pressure anf temp, pressure curve fit to gaussian 

OriginLand.py : (on left) materials; (on right) landing location 

OriginPres.py : peak pressure throughout the simulation (only ejecta plotted on right) 

OriginTemp.py : (on left) peak temperature throughout the simulation; (on right) peak temperature up to lauch point 

