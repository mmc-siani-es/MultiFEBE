Scordelis-Lo Roof validation test

Requirements: 
  - GNU/Linux environment (awk required)
  - Octave (or MATLAB, not tested)
  - Gmsh (installed)
  - multifebe (installed)

Input files:
  - README.txt
  - run_study_alt.m (Octave script for running the validation test, using symmetry planes)
  - run_study.m (Octave script for running the validation test, using B.C. on symmetry planes)
  - cylindrical_shell.geo (Gmsh script for mesh generation)
  - cylindrical_shell_alt.dat (multifebe input file, using symmetry planes)
  - cylindrical_shell.dat (multifebe input file)  
  
Execute:
> octave run_study.m

Output files:
  - *.pdf with plots of convergence results
