function sdf_labmatch,fname,ndat,labtest
;
; - - returns the number of sdf datasets in fname as ndat; prints a brief
; - - summary of header info for each dataset.
;

np=N_params()
if(np ne 3) then begin
   Message, 'Usage:',/info
   Message, 'order=sdf_labmatch,fname,ndat,labtest', /info
   Message, 'Purpose: Determine the order number of datasets with label ', /info
   Message, '         equal to labtest', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, 'Input:   labtest  = the string to test against label', /info
   Message, 'Output:  ndat  = the number of datasets in the file', /info
   Message, 'Return:  array of length ndat containing dataset orders of', /info
   Message, '         match results', /info
   Message, '         with values of -1 corresponding to no match', /info
   return, -1L
endif

sdfstr="SDF format"
get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
if(error_diag ne 0) then begin
   print, "sdf_labmatch: file probably doesn't seem to exist, returning"
   ndat=0
   close,unit
   free_lun,unit
   return, -1L
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte[0:9])

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_labmatch: fname = ",fname," is not an SDF file, returning"
   ndat=0
   close,unit
   free_lun,unit
   return, [-1L]
endif

hdrpos=0LL
datapos=0LL
norder=0L
hdrsize=0LL

readu,unit,hdrpos
readu,unit,datapos
readu,unit,norder
readu,unit,hdrsize

point_lun,-unit,curpos
header=bytarr(hdrpos-curpos)

readu,unit,header
headerstring=string(header)
;datastrings=strsplit(headerstring,string(10b),/extract)
datastrings=strtok(headerstring,string(10b),/extract)
ndat=n_elements(datastrings)
matchind=long(replicate(-1L,ndat))
ind=0L
for i=0L,long(ndat)-1L do begin
    string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
    if(labtest eq label) then begin
       matchind(ind)=long(i)
       ind=ind+1L
    endif
endfor
close,unit
free_lun,unit
return, matchind
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
