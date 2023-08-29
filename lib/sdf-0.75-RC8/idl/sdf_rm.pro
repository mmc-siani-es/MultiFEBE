pro sdf_rm,fname,quiet=quiet
;
; - - remove the sdf file fname
;
np=N_params()
if(np ne 1) then begin
   Message, 'Usage:',/info
   Message, 'sdf_rm,fname,/quiet', /info
   Message, 'Purpose: remove SDF file fname ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, 'Keyword: /quiet  - eliminate non-error output messages', /info
   return
endif
sdfstr="SDF format"
get_lun,unit
openu,unit,fname,/swap_if_little_endian,error=error_diag
if(error_diag ne 0) then begin
   if(keyword_set(quiet)) then begin
      ; do nothing
   endif else begin
      print, "sdf_rm: file doesn't seem to exist, returning"
   endelse
   close,unit
   free_lun,unit
   return
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte)

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_rm: fname = ",fname," is not an SDF file, won't remove it"
   close,unit
   free_lun,unit
   return
endif
close,unit
openu,unit,fname,/delete
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
