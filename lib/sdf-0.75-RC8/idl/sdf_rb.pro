pro sdf_rb,fname,datapos,data,elem0=elem0,nelem=nelem
np=N_params()
if(np ne 3) then begin
   Message, 'Usage:',/info
   Message, 'sdf_rb,fname,datapos,data,elem0=elem0,nelem=nelem,', /info
   Message, 'Purpose: Read data from file in large endian binary format ', /info
   Message, 'Input:   fname  = the name of the binary file to read', /info
   Message, 'Input and Output: datapos  = the starting location (bytes)', /info
   Message, '         in the file; on output, datapos = next available', /info
   Message, '         location for reading', /info
   Message, 'Output:  data  = the data to be read', /info
   Message, '         Note:  data *must* be pre-defined and dimensioned', /info
   Message, '         correctly or sdf_rb will exit without reading data', /info
   Message, '         or will read data incorrectly', /info
   Message, 'Keyword: elem0  = the initial element of data to be read', /info
   Message, '         (= 0 if not specified)', /info
   Message, 'Keyword: nelem  = the number of elements to read', /info
   Message, '         (= n_elements(data) if not specified)', /info
   Message, 'Note:    Since sdf_rb uses no metadata, you must have full', /info
   Message, '         pre-existing knowledge of the file structure to', /info
   Message, '         use sdf_rb correctly', /info
   return
endif

if(n_elements(data) eq 0LL) then begin
   print,"sdf_rb:  data is undefined, sdf_rb returning."
   return
endif

get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
if(error_diag ne 0) then begin
   print, "sdf_rb: can't open file ",fname," for reading; returning."
   close,unit
   free_lun,unit
   return
endif 

datapos64=long64(datapos)
curpos=long64(datapos64)
if(keyword_set(nelem)) then begin
   nelem64=long64(nelem)
endif else begin
   nelem64=long64(n_elements(data))
endelse

if(keyword_set(elem0)) then begin
   elem064=long64(elem0)
endif else begin
   elem064=0LL
endelse

elemlast64=long64(elem064+nelem64-1LL)
if(elemlast64 gt n_elements(data)-1) then begin
   print, "sdf_rb:  Try to write past end of array, returning"
   close,unit
   free_lun,unit
   return
endif

point_lun,unit,curpos
tempdata=data[elem064:elemlast64]
readu,unit,tempdata
data[elem064:elemlast64]=temporary(tempdata)
point_lun,-unit,curpos
datapos64=long64(curpos)
datapos=datapos64
close,unit
free_lun,unit
return
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
