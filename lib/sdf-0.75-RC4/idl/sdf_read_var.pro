function sdf_read_var,fname,labtest,quiet=quiet
np=N_params() 
if(np ne 2) then begin
   Message, 'Usage:',/info
   Message, 'data=sdf_read_var(fname,labtest,quiet=quiet)', /info
   Message, 'Purpose: Read in 1 dataset from SDF file with input label', /info
   Message, '         matching the dataset label', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         labtest = input label for desired dataset in file', /info
   Message, 'Output:  data = the returned data', /info
   Message, 'Note:    If there are repeated occurrences of datasets ', /info
   Message, '         with the same label, the *first* occurrence in ', /info
   Message, '         the file is returned ', /info
   Message, 'Note:    If keyword quiet is invoked, no non-error output', /info 
   Message, '         is printed ', /info
   return, -1
endif

match=sdf_labmatch(fname,ndat,labtest)
order=match[0]
if (order eq -1) then begin
   print,"sdf_read_var: no match for labels against ", labtest
   return, -1
endif
sdf_read,fname,order,label,data,quiet=quiet
if(n_elements(data) eq 1) then begin
   data=data[0]
endif
return, data
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
