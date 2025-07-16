pro sdf_query,fname,ndat,quiet=quiet
;
; - - returns the number of sdf datasets in fname as ndat; prints a brief
; - - summary of header info for each dataset.
;

np=N_params()
if(np ne 2) then begin
   Message, 'Usage:',/info
   Message, 'sdf_query,fname,ndat,/quiet', /info
   Message, 'Purpose: Print out summary of an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, 'Output:  ndat  = the number of datasets in the file', /info
   Message, 'Keyword: /quiet - suppress all non-error messages', /info
   return
endif

if(keyword_set(quiet)) then begin
   verbose = 0
endif else begin
   verbose = 1
endelse

sdfstr="SDF format"
get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
if(error_diag ne 0) then begin
   print, "sdf_query: file probably doesn't seem to exist, returning"
   ndat=0L
   close,unit
   free_lun,unit
   return
endif

;openr,unit,fname,/swap_if_little_endian
sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte[0:9])
;print,"sdf_query: sdftest = ",sdftest

;if (sdftest ne sdfstr) then begin
if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_query: fname = ",fname," is not an SDF file, returning"
   ndat=0
   close,unit
   free_lun,unit
   return
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
ndat=long(n_elements(datastrings))
for i=0L,ndat-1L do begin
    string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
    get_size,dims,datasize
    if(verbose eq 1) then begin 
       print,"sdf_query: iorder =",i," descriptor: ",datastrings[i],$
         " datasize = ",datasize
    endif
endfor
close,unit
free_lun,unit
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
