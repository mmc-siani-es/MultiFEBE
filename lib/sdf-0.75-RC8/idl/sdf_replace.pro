pro sdf_replace,fname,replace,ilabel,data, hinitsz=hinitsz
;
; - - replace dataset "replace" 
; - - in file fname.  
;
sdfstr="SDF format"
;hdrinit=2000LL
safety=100LL
np=N_params()
if(np ne 4) then begin
   Message, 'Usage:',/info
   Message, 'sdf_replace,fname,replace,ilabel,data, hinitsz=hinitsz', /info
   Message, 'Purpose: replace 1 dataset in an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         replace = desired sequence of dataset to replace', /info
   Message, '         ilabel = brief label for new dataset', /info
   Message, '         data = the data to replace old data', /info
   Message, '         hinitsz  -- keyword set to initial header size', /info
   Message, '         (if not set, uses default value of 2000 bytes)', /info
   return
endif

if(keyword_set(hinitsz)) then begin
   hdrinit=long64(hinitsz)
endif else begin
   hdrinit=2000LL
endelse

get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag

if(error_diag ne 0) then begin
   print, "sdf_replace: file probably doesn't exist, returning"
   close,unit
   free_lun,unit
   return
endif

sdftestbyte=bytarr(11)
readu,unit,sdftestbyte
sdftest=string(sdftestbyte)

if(not strcmp(sdftest,sdfstr,10)) then begin
   print,"sdf_replace: fname = ",fname," is not an SDF file, returning"
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
if((replace gt (norder-1L)) or (replace lt 0L)) then begin
   print,"sdf_replace: can't replace no. ",replace," in ", fname
   close,unit
   free_lun,unit
   return
endif

filepos=0LL
point_lun,-unit,filepos
filepos=long64(filepos)
header=bytarr(hdrpos-filepos)

readu,unit,header
headerstring=string(header)
datastrings=strtok(headerstring,string(10b),/extract)
ndat=n_elements(datastrings)

if(ndat ne norder) then begin
   print,"sdf_replace: error - no. strings = ",ndat," ne norder = ",norder
   close,unit
   free_lun,unit
   return
endif

if(replace gt (norder-1L)) then begin
   print,"sdf_replace: desired dataset no.",replace," too big for this file"
   close,unit
   free_lun,unit
   return
endif

if(hdrpos ge (hdrsize-safety)) then begin
   print,"sdf_replace:  header almost full, incrementing hdrsize by hdrinit"
   oldhdrsize=hdrsize
   hdrsize=hdrsize+hdrinit
   point_lun,unit,long64(strlen(sdfstr)+1+2*8+4)
   writeu,unit,hdrsize
   startloc=long64(strlen(sdfstr)+1+8+oldhdrsize)
   move_up,unit,hdrinit,startloc,datapos
   datapos=datapos+hdrinit
endif

;
point_lun,unit,filepos
curpos=long64(10+1+8+hdrsize)
;print,"sdf_replace: norder = ",norder
for i=0L,norder-1 do begin
   point_lun,unit,filepos
   if(i eq replace) then begin

      string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
      get_size,dims,datasize
      if(dtype eq 'c') then begin
         byte_old=datasize*nbpw*long64(2)
      endif else begin
         byte_old=datasize*long64(nbpw)
      endelse

      iorder=replace
      ;
      ; - - Now, get dimension information and datatypes for data:
      ;
      dsize=size(data)
      ndim=dsize[0]
      if(ndim eq 0) then begin
         dims=lon64arr(1)
         dims[0]=1
         ndim=1
      endif else begin
         dims=lon64arr(ndim)
         for ii=0L,ndim-1 do begin
            dims[ii]=dsize[ii+1]
         endfor
      endelse
      typecode=dsize[dsize[0]+1]
      if (typecode eq 4) then begin
         dtype = "f"
         nbpw=4
      endif else if (typecode eq 6) then begin
         dtype = "c"
         nbpw = 4
      endif else if ((typecode eq 3) or (typecode eq 13)) then begin
         dtype = "i"
         nbpw = 4
      endif else if (typecode eq 5) then begin
         dtype = "f"
         nbpw = 8
      endif else if (typecode eq 1) then begin
         dtype = "b"
         nbpw = 1
      endif else if ((typecode eq 2) or (typecode eq 12)) then begin
         dtype = "i"
         nbpw = 2
      endif else if (typecode eq 9) then begin
         dtype = "c"
         nbpw = 8
      endif else if ((typecode eq 14) or (typecode eq 15)) then begin
         dtype = "i"
         nbpw = 8
      endif else if(typecode eq 7) then begin
         print,"string data not allowed in sdf_replace.  Convert to byte data"
         close,unit
         free_lun,unit
         return
      endif else begin
         print, "illegal data type in sdf_replace, typecode = ",typecode
         if(typecode eq 0) then begin
            print, "Variable not defined?"
         endif
         close,unit
         free_lun,unit
         return
      endelse
      ;
      ; - - check to make sure datatype not blank
      ;
      if(dtype eq " ") then begin
         print,"sdf_replace: dtype seems blank, returning"
         close,unit
         free_lun,unit
         return
      endif
      ;
      ; - - make sure that label is blank-free:
      ;
      ilabel=strcompress(ilabel,/remove_all)
      ;
      ; - - make sure ilabel has *something* in it:
      ;
      if(strlen(ilabel) eq 0) then begin
         ilabel="NA"
      endif
      ;
      ; - - create the output string for header
      ;
      ; print,"sdf_replace: ilabel = ",ilabel
      data_to_string,outputs,iorder,ilabel,dtype,nbpw,ndim,dims
      hdrlen=strlen(outputs)
      hdrbytearr=bytarr(hdrlen+1)
      hdrbytearr[0:hdrlen-1]=byte(outputs)
      hdrbytearr[hdrlen]=10b ; add what I hope is a line-feed character at end
      ;
      ; - - write out the header string for newly replaced dataset
      ;
      point_lun,unit,filepos
      writeu,unit,hdrbytearr
      point_lun,-unit,filepos ; get new value of filepos
      filepos=long64(filepos)
      get_size,dims,datasize ; calculate size of new dataset
      if(dtype eq 'c') then begin
         byte_new=datasize*nbpw*long64(2)
      endif else begin
         byte_new=datasize*long64(nbpw)
      endelse
      ; - - now figure out if file needs to be bigger, smaller, or stay the same
      shift_up = byte_new - byte_old
      if (shift_up lt 0) then begin
         shift_dn = -shift_up
         move_dn,unit,shift_dn,curpos,datapos
      endif
      if (shift_up gt 0) then begin
         move_up,unit,shift_up,curpos,datapos
      endif
      point_lun,unit,curpos
      writeu,unit,data
      point_lun,-unit,curpos
      curpos=long64(curpos)
   endif else begin
      ; print,"sdf_replace: i = ",i," replace = ",replace
      string_to_data,datastrings[i],iorder,label,dtype,nbpw,ndim,dims
      get_size,dims,datasize
      if(dtype eq 'c') then begin
         byte_size=datasize*nbpw*long64(2)
      endif else begin
         byte_size=datasize*long64(nbpw)
      endelse
      curpos=curpos+byte_size
      ; - - write out header string for this pre-existing dataset
      data_to_string,outputs,iorder,label,dtype,nbpw,ndim,dims
      ; print,"sdf_replace:  iorder = ",iorder," outputs = ",outputs
      hdrlen=strlen(outputs)
      hdrbytearr=bytarr(hdrlen+1)
      hdrbytearr[0:hdrlen-1]=byte(outputs)
      hdrbytearr[hdrlen]=10b ; add what I hope is a line-feed character at end
      ;
      ; - - write out the header string for existing dataset
      ;
      point_lun,unit,filepos
      writeu,unit,hdrbytearr
      point_lun,-unit,filepos
      filepos=long64(filepos)
   endelse
endfor

; update values for hdrpos, datapos,norder
hdrpos=long64(filepos)
datapos=long64(curpos)
;print,"sdf_replace: new norder = ",norder
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
