// Geometry & mesh parameters
D = 1;
L = 15;
s_D = 5; 
theta = -20;  
ms_pile = 2.5;
ms_near = 1.5;
ms_far = 3;
R_truncation = 2.5*L;

theta = theta*Pi/180; 
xp = s_D/2;           
yp = s_D/2;           

Point (1) = {xp,yp,0,ms_pile};                          
Point (2) = {xp-L*Sin(theta),yp,-L*Cos(theta),ms_pile}; 
Line  (1) = {1,2};
Transfinite Line{1}=Round(L/ms_pile+1.);
Physical Line ("line-load_1") = {1};

Point (3) = {xp,yp,0,ms_pile};                           
Point (4) = {xp-L*Sin(theta),yp,-L*Cos(theta),ms_pile};   
Line  (2) = {3,4};
Transfinite Line{2}=Round(L/ms_pile+1.);
Physical Line ("pile_1") = {2};

Point (17) = {            0 ,            0 , 0 , ms_near };
Point (18) = { R_truncation ,            0 , 0 , ms_far };
Point (19) = {            0 , R_truncation , 0 , ms_far };

Line  (9) = {17,18};
Circle(10) = {18,17,19};
Line  (11) = {19,17};
Line Loop (12) = {9,10,11}; 
Plane Surface(1) = {12};

// Embedding a point in the surface 
Point (20) = {xp,yp,0,ms_near};  
Point {20} In Surface {1};

Physical Surface("free-surface") = {1};

// Mesh generation
Mesh 2;
SetOrder 2;
Save "symmetric_quarter_group.msh";
