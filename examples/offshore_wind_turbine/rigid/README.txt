To create the mesh (mesh.geo -> mesh.msh) via command-line:

$ gmsh -order 2 -2 -format msh2 mesh.geo 

To create the mesh (mesh.geo -> mesh.msh) via GUI, you have to open mesh.geo, then
mesh it (by clicking on the tree: Mesh > 2D > Set Order 2), and finally clicking
File > Export > Format: Mesh: Gmsh MSH (*.msh) > Filename: mesh.msh > Version 2 ASCII

To run the solver execute:

$ multifebe -i case.dat

Once solved, you have all the nodal results in case.dat.nso. You can visualize
in Gmsh by opening case.dat.pos. Via command-line it could be done with:

$ gmsh case.dat.pos
