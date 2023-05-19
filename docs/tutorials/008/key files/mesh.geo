// Geometry & mesh parameters
D = 10;
L = 10;
ms_near = 1;
ms_far = 2;
R_truncation = 3*D;

/*
   bucket-skirt
*/
Point (1) = { D/2 ,   0 ,  0 , ms_near };
Point (2) = {   0 ,   0 ,  0 , ms_near };
Point (3) = {   0 , D/2 ,  0 , ms_near };
Point (4) = {   0 , D/2 , -L , ms_near };
Point (5) = {   0 ,   0 , -L , ms_near };
Point (6) = { D/2 ,   0 , -L , ms_near };
Circle (1) = {1,2,3};
Line (2) = {3,4}; 
Circle (3) = {4,5,6};
Line (4) = {6,1};
Line Loop (5) = {1,2,3,4};
Transfinite Line{1,3}=Round(0.25*Pi*D/ms_near+1.);
Transfinite Line{2,4}=Round(L/ms_near+1.);
Ruled Surface (1) = {5};
Transfinite Surface {1};
Recombine Surface{1};
Color {196,196,196} { Surface {1};}

/*
   soil-skirt
*/
Point (11) = { D/2 ,   0 ,  0 , ms_near };
Point (12) = {   0 ,   0 ,  0 , ms_near };
Point (13) = {   0 , D/2 ,  0 , ms_near };
Point (14) = {   0 , D/2 , -L , ms_near };
Point (15) = {   0 ,   0 , -L , ms_near };
Point (16) = { D/2 ,   0 , -L , ms_near };
Circle (11) = {11,12,13};
Line (12) = {13,14}; 
Circle (13) = {14,15,16};
Line (14) = {16,11};
Line Loop (15) = {11,12,13,14};
Transfinite Line{11,13}=Round(0.25*Pi*D/ms_near+1.);
Transfinite Line{12,14}=Round(L/ms_near+1.);
Ruled Surface (11) = {15};
Transfinite Surface {11};
Recombine Surface{11};
Color {196,98,0} { Surface {11};}

/*
   Pile top
*/
Point (17) = {   0 ,   0 , 0 , ms_near };
Point (18) = { D/2 ,   0 , 0 , ms_near };
Point (19) = {   0 , D/2 , 0 , ms_near };
Circle(17) = {18,17,19};
Line  (18) = {19,17};
Line  (19) = {17,18};
Transfinite Line{17}=Round(0.25*Pi*D/ms_near+1.);
Line Loop (20) = {17,18,19};
Plane Surface(17) = {20};
Color {128,64,0} { Surface {17};}

/*
   Free-surface
*/
Point (31) = {   0 ,   0 , 0 , ms_near };
Point (32) = { D/2 ,   0 , 0 , ms_near };
Point (33) = {   0 , D/2 , 0 , ms_near };
Point (34) = { R_truncation ,           0. , 0. , ms_far };
Point (35) = {           0. , R_truncation , 0. , ms_far };
Circle(31) = {32,31,33};
Transfinite Line{31}=Round(0.25*Pi*D/ms_near+1.);
Circle(32) = {34,31,35};
Line  (33) = {32,34};
Line  (34) = {35,33};
Line Loop (35) = {33,32,34,-31}; 
Plane Surface(31) = {35};
Color {255,128,0} { Surface {31};}

/* Final order of Physical Surfaces */
Physical Surface("free-surface") = {31};
Physical Surface("bucket-lid") = {17};
Physical Surface("soil-skirt") = {11};
Physical Surface("bucket-skirt") = {1};

// Mesh generation
Mesh 2;
SetOrder 2;
Save "mesh.msh";
