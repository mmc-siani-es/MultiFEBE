Include "ne.geo";
If(incom == 1)
  Mesh.SecondOrderIncomplete = 1;
EndIf
ms=1;
Point(1)={0,0,0,ms};
Point(2)={1.,0,0,ms};
Point(3)={1.,1.,0,ms};
Point(4)={0,1.,0,ms};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Transfinite Line{1,2,3,4}=ne+1;
Line Loop(1)={1,2,3,4};
Surface(1)={1};
Transfinite Surface{1};
If(elem == 1)
  Recombine Surface{1};
EndIf
Physical Surface("Rectangular_plate")={1};
