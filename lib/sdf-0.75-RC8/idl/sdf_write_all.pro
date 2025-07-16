pro sdf_write_all, fname,suffix=suffix,append=append,hinitsz=hinitsz
np=N_params()
if(np ne 1) then begin
   Message, 'Usage:',/info
   Message, 'sdf_write_all,fname,suffix=suffix, /append, hinitsz=hinitsz', /info
   Message, 'Purpose: Write out all IDL vars to an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, 'Input:   valid IDL variables at calling level (level=-1)', /info
   Message, 'Output:  the SDF file, containing all the recognized IDL', /info
   Message, '         variables defined at the calling level of session', /info
   Message, 'keyword: suffix = trailing string to add to label in file', /info
   Message, 'keyword: /append = set to append data to an existing file', /info
   Message, 'keyword: hinitsz = initial header size (optional)', /info
   return
endif

if(keyword_set(append)) then begin
  ; do nothing
endif else begin
  sdf_rm,fname,/quiet
endelse

if(keyword_set(hinitsz)) then begin
  hinitstr=strcompress(',/quiet,hinitsz='+string(hinitsz),/remove_all)
endif else begin
  hinitstr=',/quiet'
endelse

var_all=strlowcase(scope_varname(level=-1))
nvar=n_elements(var_all)
newvar_all=strarr(nvar)
for i=0,nvar-1 do begin
   if(keyword_set(suffix)) then begin
      newvar_all[i]=var_all[i]+suffix
   endif else begin
      newvar_all[i]=var_all[i]
   endelse
endfor
for i=0,nvar-1 do begin
   scopestr='(scope_varfetch('+'"'+var_all[i]+'", level=-1))'+hinitstr
   ; print,'scopestr = ',scopestr
   stringcom='sdf_write,'+'"'+fname+'"'+','+'"'+newvar_all[i]+'"'+','+scopestr
   ; print,'stringcom = ',stringcom
   result=execute(stringcom)
endfor
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
