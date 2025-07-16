pro sdf_transpose,indorder,directions,data
np=N_params() 
if(np ne 3) then begin
   Message, 'Usage:',/info
   Message, 'sdf_transpose,indorder,directions,data', /info
   Message, 'Purpose: Perform multi-dimensional transpose of data', /info
   Message, '         according to the order specified by indorder', /info
   Message, 'Input:   indorder  = an array specifying the order of the',/info
   Message, '         new indices and dimensions in the transposed array', /info
   Message, 'Input:   directions  = an array specifying the direction of ',/info
   Message, '         new indices corresponding to the new dimensions in',/info 
   Message, '         the transposed array.  If any element of', /info
   Message, '         directions is < 0, the corresponding set of ', /info
   Message, '         indices is reversed', /info
   Message, 'Input:   data = the input array ', /info
   Message, 'Output:  data = the transposed array, returned in place', /info
   return
endif

arsize=size(data)
ndim=arsize[0]
if((n_elements(indorder) lt ndim) or (n_elements(directions) lt ndim)) $
then begin
   print,"sdf_transpose: not enough elements in indorder or directions"
   return
endif
if(ndim gt 1) then begin
   data=temporary(transpose(data,indorder[0:ndim-1])); use IDL transpose  
;
; - because the vacancy cycle code to do the in-place transpose is 
; - very slow in IDL
; - we use IDL's own transpose function.  Although the result is returned
; - "in place", IDL's transpose appears to create a temporary array, meaning
; - that if the input array is very large, the operation will use swap space.
;
endif
mindir=min(directions[0:ndim-1])
if(mindir ge 0) then begin
   return
endif
for i=0,ndim-1 do begin
   if(directions[i] lt 0) then begin
        data=temporary(reverse(data,i+1,/overwrite))
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
