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
Point (2) = {xp-L*Sin(theta),2.5,-L*Cos(theta),ms_pile};
Line  (1) = {1,2};
Transfinite Line{1}=Round(L/ms_pile+1.);
Physical Line ("line-load_1") = {1};

Point (3) = {xp,yp,0,ms_pile};
Point (4) = {xp-L*Sin(theta),2.5,-L*Cos(theta),ms_pile};
Line  (2) = {3,4};
Transfinite Line{2}=Round(L/ms_pile+1.);
Physical Line ("pile_1") = {2};

Point (5) = {xp,-yp,0,ms_pile};
Point (6) = {xp-L*Sin(theta),-2.5,-L*Cos(theta),ms_pile};
Line  (3) = {5,6};
Transfinite Line{3}=Round(L/ms_pile+1.);
Physical Line ("line-load_2") = {3};

Point (7) = {xp,-yp,0,ms_pile};
Point (8) = {xp-L*Sin(theta),-2.5,-L*Cos(theta),ms_pile};
Line  (4) = {7,8};
Transfinite Line{4}=Round(L/ms_pile+1.);
Physical Line ("pile_2") = {4};

Point (9) = {-xp,-yp,0,ms_pile};
Point (10) = {-xp-L*Sin(-theta),-2.5,-L*Cos(-theta),ms_pile};
Line  (5) = {9,10};
Transfinite Line{5}=Round(L/ms_pile+1.);
Physical Line ("line-load_3") = {5};

Point (11) = {-xp,-yp,0,ms_pile};
Point (12) = {-xp-L*Sin(-theta),-2.5,-L*Cos(-theta),ms_pile};
Line  (6) = {11,12};
Transfinite Line{6}=Round(L/ms_pile+1.);
Physical Line ("pile_3") = {6};

Point (13) = {-xp,yp,0,ms_pile};
Point (14) = {-xp-L*Sin(-theta),2.5,-L*Cos(-theta),ms_pile};
Line  (7) = {13,14};
Transfinite Line{7}=Round(L/ms_pile+1.);
Physical Line ("line-load_4") = {7};

Point (15) = {-xp,yp,0,ms_pile};
Point (16) = {-xp-L*Sin(-theta),2.5,-L*Cos(-theta),ms_pile};
Line  (8) = {15,16};
Transfinite Line{4}=Round(L/ms_pile+1.);
Physical Line ("pile_4") = {8};

Point (17) = {            0 ,            0 , 0 , ms_near };
Point (18) = { R_truncation ,            0 , 0 , ms_far };
Point (19) = {            0 , R_truncation , 0 , ms_far };
Point (20) = {-R_truncation ,            0 , 0 , ms_far };
Point (21) = {            0, -R_truncation , 0 , ms_far };

Line  (9) = {17,18};
Circle(10) = {18,17,19};
Line  (11) = {19,17};
Line Loop (12) = {9,10,11}; 
Plane Surface(1) = {12};

Line  (13) = {20,17};
Circle(14) = {19,17,20};
Line Loop (15) = {13,-11,14}; 
Plane Surface(2) = {15};

Line  (16) = {21,17};
Circle(17) = {20,17,21};
Line Loop (18) = {16,-13,17}; 
Plane Surface(3) = {18};

Circle(19) = {21,17,18};
Line Loop (20) = {-16,19,-9}; 
Plane Surface(4) = {20};

// Embedding points in the surfaces 
Point (22) = {xp,yp,0,ms_near};
Point (23) = {xp,-yp,0,ms_near};
Point (24) = {-xp,-yp,0,ms_near};
Point (25) = {-xp,yp,0,ms_near};
Point {22} In Surface {1};
Point {23} In Surface {2};
Point {24} In Surface {3};
Point {25} In Surface {4};

Physical Surface("free-surface") = {1,2,3,4};

// Mesh generation
Mesh 2;
SetOrder 2;
Save "whole_group.msh";