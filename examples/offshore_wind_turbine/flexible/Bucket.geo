/*
  Variables requeridas: L, D, ms_foundation, idpc,xc,yc,zc, idlid, idskirt, idbucket
*/

Function Bucket
  Geometry.AutoCoherence=0;
  // Posicion del punto intermedio en arista (0<xm<D/2)
  xm=D/4.; 
  // Coeficiente de reparto para el punto interior (provoca angulos de 120º entre lineas interiores para mejorar calidad de los elementos)
  fr=(Sqrt(3.)+1.)/(2.*Sqrt(3.)); 
  // Id central point
  poc=idpc;
  // New points
  po1 =newp; Point(po1 )={xc+xm*fr,yc+xm*fr,zc,ms_foundation};
  po2 =newp; Point(po2 )={xc-xm*fr,yc+xm*fr,zc,ms_foundation};
  po3 =newp; Point(po3 )={xc-xm*fr,yc-xm*fr,zc,ms_foundation};
  po4 =newp; Point(po4 )={xc+xm*fr,yc-xm*fr,zc,ms_foundation};
  po5 =newp; Point(po5 )={xc,yc+xm,zc,ms_foundation};
  po6 =newp; Point(po6 )={xc-xm,yc,zc,ms_foundation};
  po7 =newp; Point(po7 )={xc,yc-xm,zc,ms_foundation};
  po8 =newp; Point(po8 )={xc+xm,yc,zc,ms_foundation};
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
  Rotate {{0,0,1},{xc,yc,0},angle} { Point{po1,po2,po3,po4,po5,po6,po7,po8,po9,po10,po11,po12,po13,po14,po15,po16,po18,po19,po20,po21,po22,po23,po24,po25}; }   
  // New Lines
  li1 =newl; Line(li1 )={po1,po5};
  li2 =newl; Line(li2 )={po5,po2};
  li3 =newl; Line(li3 )={po2,po6};
  li4 =newl; Line(li4 )={po6,po3};
  li5 =newl; Line(li5 )={po3,po7};
  li6 =newl; Line(li6 )={po7,po4};
  li7 =newl; Line(li7 )={po4,po8};
  li8 =newl; Line(li8 )={po8,po1};
  li9 =newl; Line(li9 )={po8,poc};
  li10=newl; Line(li10)={po5,poc};
  li11=newl; Line(li11)={po6,poc};
  li12=newl; Line(li12)={po7,poc};
  li13=newl; Line(li13)={po9,po8};
  li14=newl; Line(li14)={po10,po5};
  li15=newl; Line(li15)={po11,po6};
  li16=newl; Line(li16)={po12,po7};
  li17=newl; Line(li17)={po13,po1};
  li18=newl; Line(li18)={po14,po2};
  li19=newl; Line(li19)={po15,po3};
  li20=newl; Line(li20)={po16,po4};
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
  Transfinite Line {li1 ,li2 ,li3 ,li4 ,li5 ,li6 ,li7 }=Ceil(0.125*Pi*D/ms_foundation)+1;
  Transfinite Line {li8 ,li9 ,li10,li11,li12,li13,li14}=Ceil(0.125*Pi*D/ms_foundation)+1;
  Transfinite Line {li15,li16,li17,li18,li19,li20,li21}=Ceil(0.125*Pi*D/ms_foundation)+1;
  Transfinite Line {li22,li23,li24,li25,li26,li27,li28}=Ceil(0.125*Pi*D/ms_foundation)+1;
  Transfinite Line {li29,li30,li31,li32,li33,li34,li35,li36}=Ceil(L/ms_foundation)+1;
  Transfinite Line {li37,li38,li39,li40,li41,li42,li43,li44}=Ceil(0.125*Pi*D/ms_foundation)+1;
  ll1 =newll; Line Loop(ll1 )={li1,li10,-li9,li8};
  ll2 =newll; Line Loop(ll2 )={li3,li11,-li10,li2};
  ll3 =newll; Line Loop(ll3 )={li5,li12,-li11,li4};
  ll4 =newll; Line Loop(ll4 )={li7,li9,-li12,li6};
  ll5 =newll; Line Loop(ll5 )={li21,li17,-li8,-li13};
  ll6 =newll; Line Loop(ll6 )={li22,li14,-li1,-li17};
  ll7 =newll; Line Loop(ll7 )={li23,li18,-li2,-li14};
  ll8 =newll; Line Loop(ll8 )={li24,li15,-li3,-li18};
  ll9 =newll; Line Loop(ll9 )={li25,li19,-li4,-li15};
  ll10=newll; Line Loop(ll10)={li26,li16,-li5,-li19};
  ll11=newll; Line Loop(ll11)={li27,li20,-li6,-li16};
  ll12=newll; Line Loop(ll12)={li28,li13,-li7,-li20};
  ll13=newll; Line Loop(ll13)={li29,li37,-li33,-li21};
  ll14=newll; Line Loop(ll14)={li33,li38,-li30,-li22};
  ll15=newll; Line Loop(ll15)={li30,li39,-li34,-li23};
  ll16=newll; Line Loop(ll16)={li34,li40,-li31,-li24};
  ll17=newll; Line Loop(ll17)={li31,li41,-li35,-li25};
  ll18=newll; Line Loop(ll18)={li35,li42,-li32,-li26};
  ll19=newll; Line Loop(ll19)={li32,li43,-li36,-li27};
  ll20=newll; Line Loop(ll20)={li36,li44,-li29,-li28};
  ss1 =news; Plane Surface(ss1 ) = {-ll1 };
  ss2 =news; Plane Surface(ss2 ) = {-ll2 };
  ss3 =news; Plane Surface(ss3 ) = {-ll3 };
  ss4 =news; Plane Surface(ss4 ) = {-ll4 };
  ss5 =news; Plane Surface(ss5 ) = {-ll5 };
  ss6 =news; Plane Surface(ss6 ) = {-ll6 };
  ss7 =news; Plane Surface(ss7 ) = {-ll7 };
  ss8 =news; Plane Surface(ss8 ) = {-ll8 };
  ss9 =news; Plane Surface(ss9 ) = {-ll9 };
  ss10=news; Plane Surface(ss10) = {-ll10};
  ss11=news; Plane Surface(ss11) = {-ll11};
  ss12=news; Plane Surface(ss12) = {-ll12};
  ss13=news; Surface(ss13) = {-ll13};
  ss14=news; Surface(ss14) = {-ll14};
  ss15=news; Surface(ss15) = {-ll15};
  ss16=news; Surface(ss16) = {-ll16};
  ss17=news; Surface(ss17) = {-ll17};
  ss18=news; Surface(ss18) = {-ll18};
  ss19=news; Surface(ss19) = {-ll19};
  ss20=news; Surface(ss20) = {-ll20};   
  Transfinite Surface {ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9,ss10,ss11,ss12};
  Recombine Surface {ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9,ss10,ss11,ss12};
  Transfinite Surface {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
  Recombine Surface {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
  Physical Surface (idlid) += {ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9,ss10,ss11,ss12};
  Physical Surface (idskirt) += {ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};
Return

