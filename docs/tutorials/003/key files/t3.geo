L = 1;  // Square side length
es = 0.15; // Element size

// Front Face
Point (1) = {0, 0, L, es};
Point (2) = {L, 0, L, es};
Point (3) = {L, L, L, es};
Point (4) = {0, L, L, es};
Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 1};
Line Loop (1) = {1, 2, 3, 4};
Plane Surface (1) = {1};
Physical Surface ("front", 1) = {1};

// Right Face
Point(5) = {L, 0, L, es};
Point(6) = {L, 0, 0, es};
Point (7) = {L, L, 0, es};
Point (8) = {L, L, L, es};
Line (5) = {5, 6};
Line (6) = {6, 7};
Line (7) = {7, 8};
Line (8) = {8, 5};
Line Loop (2) = {5, 6, 7, 8};
Plane Surface (2) = {2};
Physical Surface ("right", 2) = {2};

// Back Face
Point(9) = {L, 0, 0, es};
Point(10) = {0, 0, 0,  es};
Point (11) = {0, L, 0, es};
Point (12) = {L, L, 0, es};
Line (9) = {9, 10};
Line (10) = {10, 11};
Line (11) = {11, 12};
Line (12) = {12, 9};
Line Loop (3) = {9, 10, 11, 12};
Plane Surface (3) = {3};
Physical Surface ("back", 3) = {3};

// Left Face
Point(13) = {0, 0, 0, es};
Point(14) = {0, 0, L, es};
Point (15) = {0, L, L, es};
Point (16) = {0, L, 0, es};
Line (13) = {13, 14};
Line (14) = {14, 15};
Line (15) = {15, 16};
Line (16) = {16, 13};
Line Loop (4) = {13, 14, 15, 16};
Plane Surface (4) = {4};
Physical Surface ("left", 4) = {4};

// Top Face
Point(17) = {0, L, L, es};
Point(18) = {L, L, L, es};
Point (19) = {L, L, 0, es};
Point (20) = {0, L, 0, es};
Line (17) = {17, 18};
Line (18) = {18, 19};
Line (19) = {19, 20};
Line (20) = {20, 17};
Line Loop (5) = {17, 18, 19, 20};
Plane Surface (5) = {5};
Physical Surface ("top", 5) = {5};

// Bottom Face
Point(21) = {0, 0, L, es};
Point(22) = {0, 0, 0, es};
Point (23) = {L, 0, 0, es};
Point (24) = {L, 0, L, es};
Line (21) = {21, 22};
Line (22) = {22, 23};
Line (23) = {23, 24};
Line (24) = {24, 21};
Line Loop (6) = {21, 22, 23, 24};
Plane Surface (6) = {6};
Physical Surface ("bottom", 6) = {6};

// Mesh generation
Mesh 2;
SetOrder 2;
Save "t3.msh";
