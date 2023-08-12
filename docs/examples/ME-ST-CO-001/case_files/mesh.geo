// Geometry & mesh parameters
Lpile = 25;
ms_pile = 2;

Point (1) = {0, 0, 0, ms_pile};
Point (2) = {0, 0, Lpile, ms_pile};
Line  (1) = {1, 2};
Transfinite Line {1} = Ceil(Lpile/ms_pile + 1.);
Physical Line ("line-load") = {1};

Point (3) = {0, 0, 0, ms_pile};
Point (4) = {0, 0, Lpile, ms_pile};
Line  (2) = {3, 4};
Transfinite Line {2} = Ceil(Lpile/ms_pile + 1.);
Physical Line ("pile") = {2};

// Mesh generation
Mesh 1;
SetOrder 2;
Save "mesh.msh";