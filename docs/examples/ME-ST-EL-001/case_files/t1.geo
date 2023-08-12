L = 1;  // Square side length
es = 1; // Element size

// Bottom
Point (1) = {0, 0, 0, es};
Point (2) = {L, 0, 0, es};
Line (1) = {1, 2};
Physical Line ("bottom", 1) = {1};

// Right
Point (3) = {L, 0, 0, es};
Point (4) = {L, L, 0, es};
Line (2) = {3, 4};
Physical Line ("right", 2) = {2};

// Top
Point (5) = {L, L, 0, es};
Point (6) = {0, L, 0, es};
Line (3) = {5, 6};
Physical Line ("top", 3) = {3};

// Left
Point (7) = {0, L, 0, es};
Point (8) = {0, 0, 0, es};
Line (4) = {7, 8};
Physical Line ("left", 4) = {4};

// Mesh generation
Mesh 1;
SetOrder 1;
Save "t1.msh";