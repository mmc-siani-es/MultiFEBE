pro sdf_wb,fname,datapos,data,elem0=elem0,nelem=nelem
np=N_params()
if(np ne 3) then begin
   Message, 'Usage:',/info
   Message, 'sdf_wb,fname,datapos,data,elem0=elem0,nelem=nelem', /info
   Message, 'Purpose: Write data to file in large endian binary format ', /info
   Message, 'Input:   fname  = the name of the file for binary I/O', /info
   Message, 'Input and Output: datapos  = the location (bytes) to start', /info
   Message, '         writing in the file on input; on output, datapos ', /info
   Message, '         = next available location for writing.', /info
   Message, '         Be sure to set datapos to 0LL for initial call', /info
   Message, 'Input:   data  = the data to be written', /info
   Message, '         The array data must intially exist', /info
   Message, 'Keyword: elem0  = the initial element of data to be written', /info
   Message, '         (= 0 if not specified)', /info
   Message, 'Keyword: nelem  = the number of elements to write', /info
   Message, '         (= n_elements(data) if not specified)', /info
   return
endif

if(n_elements(data) eq 0LL) then begin
   print,"sdf_wb:  data is not defined, returning"
   return
endif

get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
new_file=0
if(error_diag ne 0) then begin
   new_file=1
   close,unit
endif else begin
   close,unit
endelse

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
   print, "sdf_wb:  Try to write past end of array, returning"
   free_lun,unit
   return
endif

if(new_file eq 0) then begin
   openu,unit,fname,/swap_if_little_endian,error=error_diag
endif else begin
   openw,unit,fname,/swap_if_little_endian,error=error_diag
endelse

if(error_diag ne 0) then begin
   print, "sdf_wb:  Problem opening file ",fname," for writing, returning"
   close,unit
   free_lun,unit
   return
endif

point_lun,unit,curpos
writeu,unit,data[elem064:elemlast64]
point_lun,-unit,curpos
datapos64=long64(curpos)
datapos=datapos64
truncate_lun,unit
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
