This folder contains an example of a seismic analysis of an offshore wind turbine.

The folder "flexible" contains the BEM-FEM model considering direct modelling of 
soil-structure, where suction caissons and soil are modelled via FEM and BEM. 
Seismic input is a vertically incident SH wave with unitary displacements
at the free-surface level. 

The folder "rigid" contains the FEM model of the offshore wind turbine. In this
case, no soil-structure interaction is considered, and unitary horizontal 
displacements are directly applied at the bottom leg nodes.

The file "plot.gp" contains a GNUPLot script for plotting the horizontal 
displacement of the RNA.
