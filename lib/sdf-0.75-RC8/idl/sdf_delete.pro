pro sdf_delete,fname,idelete
;
; - - deletes dataset idelete 
; - - from file fname.  
;
sdfstr="SDF format"
;hdrinit=2000LL
safety=100LL
np=N_params()
if(np ne 2) then begin
   Message, 'Usage:',/info
   Message, 'sdf_delete,fname,idelete', /info
   Message, 'Purpose: Delete 1 dataset of an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         idelete = sequence of dataset to delete', /info
   return
endif
get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag

if(error_diag ne 0) then begin
   print, "sdf_delete: file probably doesn't exist, returning"
   close,unit
   free_lun,unit
   return
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte)

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_delete: fname = ",fname," is not an SDF file, returning"
   close,unit
   free_lun,unit
   return
endif
close,unit

hdrpos=0LL
datapos=0LL
norder=0L
hdrsize=0LL


; - get current values of header position, data postion, and dataset order:

openu,unit,fname,/swap_if_little_endian
point_lun,unit,long64(11)
readu,unit,hdrpos
readu,unit,datapos ; datapos is also the current filesize
readu,unit,norder
readu,unit,hdrsize
if((idelete ge norder) or (idelete lt 0L)) then begin
   print,"sdf_delete: can't delete non-existent dataset no. ",idelete," in ",$
	   fname
   close,unit
   free_lun,unit
   return
endif
if(norder eq 1) then begin
   close,unit
   free_lun,unit
   sdf_rm,fname
   return
endif

filepos=0LL
point_lun,-unit,filepos
filepos=long64(filepos)
header=bytarr(hdrpos-filepos)

readu,unit,header
headerstring=string(header)
;datastrings=strsplit(headerstring,string(10b),/extract)
datastrings=strtok(headerstring,string(10b),/extract)
ndat=n_elements(datastrings)

if(ndat ne norder) then begin
   print,"sdf_delete: error - no. strings = ",ndat," ne norder = ",norder
   close,unit
   free_lun,unit
   return
endif

if(idelete ge norder) then begin
   print,"sdf_delete: desired dataset no.",idelete," too big for this file"
   close,unit
   free_lun,unit
   return
endif
;
;filepos=0LL
;point_lun,-unit,filepos
;header=bytarr(hdrpos-filepos)
;
;readu,unit,header
;headerstring=string(header)
;datastrings=strsplit(headerstring,string(10b),/extract)
;ndat=n_elements(datastrings)
;
curpos=long64(10+1+8+hdrsize)
;print,"sdf_delete: norder = ",norder
for i=0L,norder-1L do begin
   string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
   if(i ge idelete) then begin
;      iorder=i-2L ; why 2???
       iorder=iorder-1
   endif
   get_size,dims,datasize
;  print,"sdf_delete: iorder = ",iorder
;  datastrings[i]=""
;     print,"sdf_delete-1: datastrings[",i,"] = ",datastrings[i]
   data_to_string,tempstring,iorder,label,dtype,nbpw,ndim,dims
   datastrings[i]=tempstring
;     print,"sdf_delete0: datastrings[",i,"] = ",datastrings[i]
   if(dtype eq 'c') then begin
      byte_size=datasize*nbpw*long64(2)
   endif else begin
      byte_size=datasize*long64(nbpw)
   endelse
   if (i eq idelete) then begin ; this is the dataset to be deleted
      move_dn,unit,byte_size,curpos,datapos
      point_lun,unit,filepos
   endif else begin
      hdrlen=strlen(datastrings[i])
      hdrbytearr=bytarr(hdrlen+1)
      hdrbytearr[0:hdrlen-1]=byte(datastrings[i])
      hdrbytearr[hdrlen]=10b ; add what I hope is a line-feed character at end
      point_lun,unit,filepos
;     print,"sdf_delete1: datastrings[",i,"] = ",datastrings[i]
      writeu,unit,hdrbytearr
      point_lun,-unit,filepos ; get new value of filepos
      filepos=long64(filepos)
      curpos=curpos+byte_size
   endelse
endfor
; update values for hdrpos, datapos,norder
hdrpos=long64(filepos)
datapos=long64(curpos)
norder=norder-1
; output the new values of hdrpos,datapos,norder to the file
point_lun,unit,long64(11)
writeu,unit,long64(hdrpos)
writeu,unit,long64(datapos)
writeu,unit,long(norder)
writeu,unit,long64(hdrsize)
; - close up and return
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
