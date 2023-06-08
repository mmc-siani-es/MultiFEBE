
*IntFormat "%12i"
*RealFormat "%25.16e"
[parts]
*Set var nLayer=0
*loop layers
*Set var nLayer=LayerNum
*end layers
*nLayer
*loop layers
*LayerNum *LayerName
*end layers
[nodes]
*nPoin
*set elems(all)
*loop nodes
*NodesNum *NodesCoord
*end nodes
[elements]
*set Elems(All)
*nElem
*set Elems(Linear)
*loop elems
*if(ElemsNnode==2)
*ElemsNum line2 1 *ElemsLayerNum *ElemsConec
*endif
*if(ElemsNnode==3)
*ElemsNum line3 1 *ElemsLayerNum *ElemsConec
*endif
*end elems
*set Elems(Triangle)
*loop elems
*if(ElemsNnode==3)
*ElemsNum tri3 1 *ElemsLayerNum *ElemsConec
*endif
*if(ElemsNnode==6)
*ElemsNum tri6 1 *ElemsLayerNum *ElemsConec
*endif
*end elems
*set Elems(Quadrilateral)
*loop elems
*if(ElemsNnode==4)
*ElemsNum quad4 1 *ElemsLayerNum *ElemsConec
*endif
*if(ElemsNnode==8)
*ElemsNum quad8 1 *ElemsLayerNum *ElemsConec
*endif
*if(ElemsNnode==9)
*ElemsNum quad9 1 *ElemsLayerNum *ElemsConec
*endif
*end elems


