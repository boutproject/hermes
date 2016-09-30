ISTTOK turbulence simulation
============================

To generate the input grids, run the "circle.py" script:

 $ python circle.py

The input file is for a q=7, 5eV edge temperature. This 
should be used for the initial simulation transients.
The particle and heat sources are adaptive, and use a 
PI controller to attempt to match the input density and 
temperature profiles.

To increase the temperature, restart using the 10eV then 20eV
input grids.

 
