pro data_to_string, outputs,io,label,datatype,nbpw,ndim,dims
;
; - - put the data into its string descriptor
;
ios=strcompress(string(io),/remove_all)
labels=strcompress(string(label),/remove_all)
dts=strcompress(string(datatype),/remove_all)
nbpws=strcompress(string(nbpw),/remove_all)
ndims=strcompress(string(ndim),/remove_all)
dimss=sindgen(ndim)
for i=0,ndim-1 do begin
   dimss[i]=strcompress(string(dims[ndim-1-i]),/remove_all)
endfor
outputs=""
blank=" "
outputs=ios+blank+labels+blank+dts+blank+nbpws+blank+ndims
for i=0,ndim-1 do begin
  outputs=outputs+blank+dimss[i] ; flip order for IDL
endfor
end
; Simple Data Format
; http://solarmuri.ssl.berkeley.edu/~fisher/public/software/SDF
; Copyright (C) 2006 University of California
;
; This is free software; you can redistribute it and/or
; modify it under the terms of the GNU Lesser General Public
; License as published by the Free Software Foundation;
; either version 2.1 of the License, or (at your option) any later version.
;
; This software is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
; See the GNU Lesser General Public License for more details.
;
; To view the GNU Lesser General Public License visit
; http://www.gnu.org/copyleft/lesser.html
; or write to the Free Software Foundation, Inc.,
; 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
