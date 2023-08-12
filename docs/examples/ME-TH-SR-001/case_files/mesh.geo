// Geometry & mesh parameters
Lbeam = 10;
ms_beam = 2.5;

Point (1) = {0, 0, 0, ms_beam};
Point (2) = {Lbeam, 0, 0, ms_beam};
Line  (1) = {1, 2};
Physical Line ("beam") = {1};

// Mesh generation
Mesh 1;
SetOrder 2;
Save "mesh.msh";