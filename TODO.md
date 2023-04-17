# MultiFEBE

### Todo

- New output files with results split by type of variable (displacements, tractions, nodal forces, etc...), and by physical entity, be boundary, be bodyload, be subregion, etc...
- Consider text output files in csv format.
- Since discontinuous boundary elements are partially implemented, remove any of the functionalities for these
- Implement proper point loads for acoustics and elastodynamics in BEM regions
- Change input file format for boundary conditions to a more readable one
- Implement read/write to GiD 
- Implement read/write to MSH 4 file format of Gmsh.
- Implement 3D finite elements
- Put into functions the properties transformation between elastic and poroelastic regions/materials/models
- Generalize free-term driver routine for general point of collocation


### In Progress ···

- Implementing rigid regions with coupling with BE and FE nodes
- Use type material(*)%* for everything instead of explicit properties in region data structure

