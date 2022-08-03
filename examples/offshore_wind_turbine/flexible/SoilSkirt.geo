/*
  Variables requeridas: L, D, ms_foundation, xc, yc, zc, idbodyload
*/

Function SoilSkirt
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
  po17=newp; Point(po17)={xc,yc,zc+L,ms_foundation};
  po18=newp; Point(po18)={xc+D/2.,yc,zc+L,ms_foundation};
  po19=newp; Point(po19)={xc,yc+D/2.,zc+L,ms_foundation};
  po20=newp; Point(po20)={xc-D/2.,yc,zc+L,ms_foundation};
  po21=newp; Point(po21)={xc,yc-D/2.,zc+L,ms_foundation};
  po22=newp; Point(po22)={xc+D/2.*Sqrt(2.)/2.,yc+D/2.*Sqrt(2.)/2.,zc+L,ms_foundation};
  po23=newp; Point(po23)={xc-D/2.*Sqrt(2.)/2.,yc+D/2.*Sqrt(2.)/2.,zc+L,ms_foundation};
  po24=newp; Point(po24)={xc-D/2.*Sqrt(2.)/2.,yc-D/2.*Sqrt(2.)/2.,zc+L,ms_foundation};
  po25=newp; Point(po25)={xc+D/2.*Sqrt(2.)/2.,yc-D/2.*Sqrt(2.)/2.,zc+L,ms_foundation};
  // Rotate points to orientate the lines to the center
  angle=Atan2(yc,xc);  
  Rotate {{0,0,1},{xc,yc,0},angle} { Point{po9,po10,po11,po12,po13,po14,po15,po16,po17,po18,po19,po20,po21,po22,po23,po24,po25}; }     
  // New Lines
  li21=newl; Circle(li21)={po9,poc,po13};
  li22=newl; Circle(li22)={po13,poc,po10};
  li23=newl; Circle(li23)={po10,poc,po14};
  li24=newl; Circle(li24)={po14,poc,po11};
  li25=newl; Circle(li25)={po11,poc,po15};
  li26=newl; Circle(li26)={po15,poc,po12};
  li27=newl; Circle(li27)={po12,poc,po16};
  li28=newl; Circle(li28)={po16,poc,po9};
  li29=newl; Line(li29)={po9 ,po18};
  li30=newl; Line(li30)={po10,po19};
  li31=newl; Line(li31)={po11,po20};
  li32=newl; Line(li32)={po12,po21};
  li33=newl; Line(li33)={po13,po22};
  li34=newl; Line(li34)={po14,po23};
  li35=newl; Line(li35)={po15,po24};
  li36=newl; Line(li36)={po16,po25};
  li37=newl; Circle(li37)={po18,po17,po22};
  li38=newl; Circle(li38)={po22,po17,po19};
  li39=newl; Circle(li39)={po19,po17,po23};
  li40=newl; Circle(li40)={po23,po17,po20};
  li41=newl; Circle(li41)={po20,po17,po24};
  li42=newl; Circle(li42)={po24,po17,po21};
  li43=newl; Circle(li43)={po21,po17,po25};
  li44=newl; Circle(li44)={po25,po17,po18};
  Transfinite Line {li21,li22,li23,li24,li25,li26,li27,li28}=Ceil(0.125*Pi*D/ms_foundation)+1;
  Transfinite Line {li29,li30,li31,li32,li33,li34,li35,li36}=Ceil(L/ms_foundation)+1;
  Transfinite Line {li37,li38,li39,li40,li41,li42,li43,li44}=Ceil(0.125*Pi*D/ms_foundation)+1;
  ll13=newll; Line Loop(ll13)={li29,li37,-li33,-li21};
  ll14=newll; Line Loop(ll14)={li33,li38,-li30,-li22};
  ll15=newll; Line Loop(ll15)={li30,li39,-li34,-li23};
  ll16=newll; Line Loop(ll16)={li34,li40,-li31,-li24};
  ll17=newll; Line Loop(ll17)={li31,li41,-li35,-li25};
  ll18=newll; Line Loop(ll18)={li35,li42,-li32,-li26};
  ll19=newll; Line Loop(ll19)={li32,li43,-li36,-li27};
  ll20=newll; Line Loop(ll20)={li36,li44,-li29,-li28};
  ss13=news; Surface(ss13) = {-ll13};
  ss14=news; Surface(ss14) = {-ll14};
  ss15=news; Surface(ss15) = {-ll15};
  ss16=news; Surface(ss16) = {-ll16};
  ss17=news; Surface(ss17) = {-ll17};
  ss18=news; Surface(ss18) = {-ll18};
  ss19=news; Surface(ss19) = {-ll19};
  ss20=news; Surface(ss20) = {-ll20};   
  Transfinite Surface {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
  Recombine Surface {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
  Physical Surface (idbodyload) += {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
Return

