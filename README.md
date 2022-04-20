# Astrolib

Astrolib is an astrodynamics library containing methods for vector math, coordinate system conversions, time system conversions, and the Gauss Initial Orbit Determination algorithm implementing both Gibbs and Herrick-Gibbs methods to obtain velocity. The Gauss IOD approximates an orbit using angles-only measurements, so the intent is to provide amateur astronomers with a quick way to compute satellite orbits given their optical observations. Unit testing is still needed for the main IOD method, and additional error handling is required.

A sample program shows how to run the IOD by inputting a text file with observation data. The observation data includes:

date lat lon alt az el

The date must be in unix time (as this is how time is displayed with measurements from online data centers) however future iterations will use date time, or julian date.
