Scordelis-Lo Roof validation test

Requirements: 
  - GNU/Linux environment (awk required)
  - Octave (or MATLAB, not tested)
  - Gmsh (installed)
  - multifebe (installed)

Input files:
  - README.txt

  - run_study.m (Octave script for running the validation test)
  - Rectangular_plate.geo (Gmsh script for mesh generation)
  - Rectangular_plate.dat (multifebe input file)  
  
Execute:
> octave run_study.m

Output files:
  - *.pdf with plots of convergence results
