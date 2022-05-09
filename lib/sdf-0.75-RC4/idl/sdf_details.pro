pro sdf_details,fname,dorder,label,dtype,nbpw,ndim,dims,quiet=quiet
;
; - - reads in the header info for dataset of order "dorder", and returns the
; - - string label into the label calling argument.  Also returns details
; - - about the dataset
;
sdfstr="SDF format"
;hdrinit=2000LL
safety=100LL
np=N_params()
if(np ne 7) then begin
   Message, 'Usage:',/info
   Message, 'sdf_details,fname,dorder,label,dtype,nbpw,ndim,dims', /info
   Message, 'Purpose: Read in info about 1 dataset from an SFD file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         dorder = order of desired dataset in file', /info
   Message, 'Output:  label = brief string descriptor of dataset', /info
   Message, 'Output:  dtype = type of data stored in dataset', /info
   Message, 'Output:  nbpw = No. bytes per word in dataset', /info
   Message, 'Output:  ndim = No. dimensions of dataset', /info
   Message, 'Output:  dims = array[ndim] containing data dimensions', /info
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
   print, "sdf_details: file probably doesn't seem to exist, returning"
   ndat=0
   close,unit
   free_lun,unit
   return
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte)
;print,"sdf_details: sdftest = ",sdftest

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_details: fname = ",fname," is not an SDF file, exiting"
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
;print,"sdf_details: hdrsize = ",hdrsize

point_lun,-unit,curpos
curpos=long64(curpos)
header=bytarr(hdrpos-curpos)

readu,unit,header
headerstring=string(header)
;datastrings=strsplit(headerstring,string(10b),/extract)
datastrings=strtok(headerstring,string(10b),/extract)
ndat=long(n_elements(datastrings))

if(ndat ne norder) then begin
   print,"sdf_details: error - no. strings ne norder"
   close,unit
   free_lun,unit
   return
endif
if(dorder ge norder) then begin
   print,"sdf_details: desired dataset no.",dorder," too big for this file"
   close,unit
   free_lun,unit
   return
endif
if(dorder lt 0L) then begin
   print,"sdf_details: desired dataset no.",dorder," is less than 0"
   close,unit
   free_lun,unit
   return
endif

curpos=0LL
curpos=long64(10+1+8+hdrsize)
for i=0L,dorder do begin
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
   print,"sdf_details: dataset ",iorder," descriptor ",datastrings[iorder], $
      " datasize = ", datasize
endif

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
