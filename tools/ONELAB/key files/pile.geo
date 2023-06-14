
Include "pile_common.pro";

// Geometry & mesh parameters
D = Diameter;
L = PL;
theta = TL;
ms_pile = PM;
ms_near = NM;
ms_far = FM;
R_truncation = Radius;

theta = TL*Pi/180;

Point (1) = {0,0,0,ms_pile};
Point (2) = {-L*Sin(theta),0,-L*Cos(theta),ms_pile};
Line  (1) = {1,2};
Transfinite Line{1}=Round(L/ms_pile+1.);
Physical Line ("line-load") = {1};

Point (3) = {0,0,0,ms_pile};
Point (4) = {-L*Sin(theta),0,-L*Cos(theta),ms_pile};
Line  (2) = {3,4};
Transfinite Line{2}=Round(L/ms_pile+1.);
Physical Line ("pile") = {2};

Point (5) = {            0 ,            0 , 0 , ms_near };
Point (6) = { R_truncation ,            0 , 0 , ms_far };
Point (7) = {            0 , R_truncation , 0 , ms_far };
Point (8) = {-R_truncation ,            0 , 0 , ms_far };
Point (9) = {            0, -R_truncation , 0 , ms_far };

Line  (3) = {5,6};
Circle(4) = {6,5,7};
Line  (5) = {7,5};
Line Loop (6) = {3,4,5}; 
Plane Surface(1) = {6};

Line  (7) = {8,5};
Circle(8) = {7,5,8};
Line Loop (9) = {7,-5,8}; 
Plane Surface(2) = {9};

Line  (10) = {9,5};
Circle(11) = {8,5,9};
Line Loop (12) = {10,-7,11}; 
Plane Surface(3) = {12};

Circle(13) = {9,5,6};
Line Loop (14) = {-10,13,-3}; 
Plane Surface(4) = {14};

Physical Surface("free-surface") = {1,2,3,4};

// Mesh generation
Mesh 2;
SetOrder 2;
Save "pile.msh";