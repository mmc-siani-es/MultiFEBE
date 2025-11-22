Include "ne.geo";
If(incom == 1)
  Mesh.SecondOrderIncomplete = 1;
EndIf
ms=1;
Point(1)={0,0,0,ms};
Point(2)={25.*Sin(40.*Pi/180.),0,25.*Cos(40.*Pi/180.),ms};
Point(3)={0,0,25.,ms};
Point(4)={0,25.,0,ms};
Point(5)={25.*Sin(40.*Pi/180.),25.,25.*Cos(40.*Pi/180.),ms};
Point(6)={0,25.,25.,ms};
Line  (1)={2,5};
Circle(2)={5,4,6};
Line  (3)={6,3};
Circle(4)={3,1,2};
Transfinite Line{1,2,3,4}=ne+1;
Line Loop(1)={1,2,3,4};
Surface(1)={1};
Transfinite Surface{1};
If(elem == 1)
  Recombine Surface{1};
EndIf
Physical Surface("cylindrical_shell_roof")={1};
