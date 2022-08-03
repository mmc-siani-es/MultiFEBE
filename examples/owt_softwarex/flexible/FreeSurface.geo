/*
  Variables requeridas: D, ms_foundation, xc, yc, zc
*/
Function FreeSurface
  Geometry.AutoCoherence=0;
  // Posicion del punto intermedio en arista (0<xm<D/2)
  xm=D/4.; 
  // Coeficiente de reparto para el punto interior (provoca angulos de 120º entre lineas interiores para mejorar calidad de los elementos)
  fr=(Sqrt(3.)+1.)/(2.*Sqrt(3.)); 
  // New points
  poc =newp; Point(poc )={xc,yc,zc,ms_foundation};
  
  po9 =newp; Point(po9 )={xc+D/2.,yc,zc,ms_foundation};
  po10=newp; Point(po10)={xc,yc+D/2.,zc,ms_foundation};
  po11=newp; Point(po11)={xc-D/2.,yc,zc,ms_foundation};
  po12=newp; Point(po12)={xc,yc-D/2.,zc,ms_foundation};
  
  po13=newp; Point(po13)={xc+D/2.*Sqrt(2.)/2.,yc+D/2.*Sqrt(2.)/2.,zc,ms_foundation};
  po14=newp; Point(po14)={xc-D/2.*Sqrt(2.)/2.,yc+D/2.*Sqrt(2.)/2.,zc,ms_foundation};
  po15=newp; Point(po15)={xc-D/2.*Sqrt(2.)/2.,yc-D/2.*Sqrt(2.)/2.,zc,ms_foundation};
  po16=newp; Point(po16)={xc+D/2.*Sqrt(2.)/2.,yc-D/2.*Sqrt(2.)/2.,zc,ms_foundation};
  // Rotate points to orientate the lines to the center
  angle=Atan2(yc,xc);  
  Rotate {{0,0,1},{xc,yc,0},angle} { Point{po9,po10,po11,po12,po13,po14,po15,po16}; }     
  // New Lines
  li21=newl; Circle(li21)={po9,poc,po13};
  li22=newl; Circle(li22)={po13,poc,po10};
  li23=newl; Circle(li23)={po10,poc,po14};
  li24=newl; Circle(li24)={po14,poc,po11};
  li25=newl; Circle(li25)={po11,poc,po15};
  li26=newl; Circle(li26)={po15,poc,po12};
  li27=newl; Circle(li27)={po12,poc,po16};
  li28=newl; Circle(li28)={po16,poc,po9};
  Transfinite Line {li21,li22,li23,li24,li25,li26,li27,li28}=Ceil(0.125*Pi*D/ms_foundation)+1;
  
  ll1 =newll; Line Loop(ll1 )={li21,li22,li23,li24,li25,li26,li27,li28};

Return

