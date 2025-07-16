pro string_to_data,inputs,io,label,datatype,nbpw,ndim,dims
;
; - - break apart 'inputs' string into its components
;
;tokens=strsplit(inputs,/extract)
tokens=strtok(inputs,/extract)
nt=n_elements(tokens)
ii=0
for i=0,nt-1 do begin
p=tokens[i]
if (i eq 0) then io=long(p)
if (i eq 1) then label=string(p)
if (i eq 2) then datatype=string(p)
if (i eq 3) then nbpw=long(p)
if (i eq 4) then begin
   ndim=long(p)
   dims=lon64arr(ndim)
endif
if (i ge 5) then begin
   dims[ndim-1-ii]=long64(p)
   ii=ii+1
endif
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
