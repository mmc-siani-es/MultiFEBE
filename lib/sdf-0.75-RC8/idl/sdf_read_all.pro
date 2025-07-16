pro sdf_read_all, fname,suffix=suffix,quiet=quiet
np=N_params()
if(np ne 1) then begin
   Message, 'Usage:',/info
   Message, 'sdf_read_all,fname,suffix=suffix,quiet=quiet', /info
   Message, 'Purpose: Read in all datasets from an SDF file ', /info
   Message, '         assuming variable names = label', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, 'Output:  all of the variables from the SDF file', /info
   Message, '         are placed into level -1 (calling level) of IDL', /info
   Message, 'Note:    if there are repeated occurrences of a dataset ', /info
   Message, '         with the same label, the *last* occurrence ', /info
   Message, '         of the dataset with a given label is returned ', /info
   Message, 'keyword: suffix - additional trailing string for var names', /info
   Message, 'keyword: keyword quiet suppresses non-error output ', /info
   return
endif

nvar=0L
sdf_query,fname,nvar,quiet=quiet
for i=0L,long(nvar)-1L do begin
   stri=strcompress(string(i))
   sdf_details,fname,i,label,dtype,nbpw,ndim,dims,quiet=quiet
   if(keyword_set(suffix)) then begin
      newlabel=label+suffix
   endif else begin
      newlabel=label
   endelse
   scopestr='(scope_varfetch('+'"'+newlabel+'", /enter, level=-1))'
   stringcom='sdf_read,'+'"'+fname+'"'+','+stri+',labelagain,'$
     +scopestr+',/quiet'
   result=execute(stringcom)
   scopestr='(scope_varfetch('+'"'+newlabel+'", level=-1))'
   oneelement='if(n_elements('+scopestr+') eq 1) then '+scopestr+$
       '='+scopestr+'[0]'
   result=execute(oneelement)
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
