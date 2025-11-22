Spherical dome validation test

Requirements: 
  - GNU/Linux environment (awk required)
  - Octave (or MATLAB, not tested)
  - Gmsh (installed)
  - multifebe (installed)

Input files:
  - README.txt
  - run_study.m (Octave script for running the validation test)
  - spherical_dome.geo (Gmsh script for mesh generation)
  - spherical_dome.dat (multifebe input file)
  
Execute:
> octave run_study.m

Output files:
  - *.pdf with plots of convergence results
