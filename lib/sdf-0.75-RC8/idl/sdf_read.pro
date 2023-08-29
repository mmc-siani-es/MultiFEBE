pro sdf_read,fname,dorder,label,data,quiet=quiet
;
; - - reads in the dataset of order "dorder", and returns the
; - - data into the data array.
; - - summary of header info for each dataset.
;
sdfstr="SDF format"
;hdrinit=2000LL
safety=100LL
np=N_params()
if(np ne 4) then begin
   Message, 'Usage:',/info
   Message, 'sdf_read,fname,dorder,label,data, /quiet', /info
   Message, 'Purpose: Read in 1 dataset from an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         dorder = order of desired dataset in file', /info
   Message, 'Output:  label = brief string descriptor of dataset', /info
   Message, '         data  = the returned data', /info
   Message, '         /quiet - set to suppress printing non-errors msgs ', /info
   return
endif

if(keyword_set(quiet)) then begin
   verbose = 0
endif else begin
   verbose = 1
endelse

get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
if(error_diag ne 0) then begin
   print, "sdf_read: file probably doesn't exist, returning"
   close,unit
   free_lun,unit
   return
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte)
;print,"sdf_read: sdftest = ",sdftest

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_read: fname = ",fname," is not an SDF file, exiting"
   close,unit
   free_lun,unit
   return
endif

hdrpos=0LL
datapos=0LL
norder=0L
hdrsize=0LL

readu,unit,hdrpos
readu,unit,datapos ; datapos is also the current filesize
readu,unit,norder
readu,unit,hdrsize
;print,"sdf_read: hdrsize = ",hdrsize

point_lun,-unit,curpos
curpos=long64(curpos)
header=bytarr(hdrpos-curpos)

readu,unit,header
headerstring=string(header)
;datastrings=strsplit(headerstring,string(10b),/extract)
datastrings=strtok(headerstring,string(10b),/extract)
ndat=n_elements(datastrings)

if(ndat ne norder) then begin
   print,"sdf_read: error - no. strings ne norder"
   close,unit
   free_lun,unit
   return
endif
if(dorder ge norder) then begin
   print,"sdf_read: desired dataset no.",dorder," too big for this file"
   close,unit
   free_lun,unit
   return
endif
if(dorder lt 0L) then begin
   print,"sdf_read: desired dataset no.",dorder," is less than 0"
   close,unit
   free_lun,unit
   return
endif

curpos=0LL
curpos=long64(10+1+8+hdrsize)
for i=0,dorder do begin
    string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
    get_size,dims,datasize
    if(dtype ne "c") then begin ; test for complex var
	nbytes=datasize*long64(nbpw)
    endif else begin
        nbytes=2LL*datasize*long64(nbpw) ; kludge for complex var
    endelse
    if(i lt dorder) then curpos=curpos+nbytes
endfor

if(verbose eq 1) then begin
    print,"sdf_read: dataset ",iorder," descriptor ",datastrings[iorder], $
      " datasize = ", datasize
endif

if((dtype eq "f") and (nbpw eq 4)) then begin
    data=fltarr(dims)
endif else if ((dtype eq "f") and (nbpw eq 8)) then begin
    data=dblarr(dims)
endif else if ((dtype eq "i") and (nbpw eq 4)) then begin
    data=lonarr(dims)
endif else if ((dtype eq "i") and (nbpw eq 8)) then begin
    data=lon64arr(dims)
endif else if ((dtype eq "c") and (nbpw eq 4)) then begin
    data=complexarr(dims)
endif else if ((dtype eq "c") and (nbpw eq 8)) then begin
    data=dcomplexarr(dims)
endif else if ((dtype eq "i") and (nbpw eq 2)) then begin
    data=intarr(dims)
endif else if((nbpw eq 1)) then begin
    data=bytarr(dims)
endif else begin
    print,"sdf_read:  Unrecognized data type: cannot create data array"
    close,unit
    free_lun,unit
    return
endelse

point_lun,unit,curpos
readu,unit,data
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
