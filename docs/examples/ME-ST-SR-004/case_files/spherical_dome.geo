Include "ne.geo";
If(incom == 1)
  Mesh.SecondOrderIncomplete = 1;
EndIf

// Data
R=10;             // Sphere radius
beta=18.*Pi/180.; // Top hole angle

// Arc at y=0
Point(1)={0,0,0,1.};
Point(2)={R,0,0,1.};
Point(3)={R*Sin(beta),0,R*Cos(beta),1.};
Circle(1)={2,1,3};

// Extrusion around z axis
out[] = Extrude { {0,0,1}, {0,0,0}, Pi/2. } { Line{1}; };
// out[0]: final line of the extrusion
// out[1]: surface
// out[2]: side line 1
// out[3]: side line 2 

// Transfinite mesh
Transfinite Line{1,out[0],out[2],out[3]}=ne+1;
Transfinite Surface{out[1]};

// Type of element
If(elem == 1)
  Recombine Surface{out[1]};
EndIf

// Physical part
Physical Surface("spherical_dome")={out[1]};

