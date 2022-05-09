pro sdf_write,fname,label,data,quiet=quiet,hinitsz=hinitsz
;
; - - writes one dataset contained in the input
; - - array data into the file fname.  The description of the data is
; - - contained within the label fields.
; - - which the user inputs.  The datatype, nbpw, and dimension information 
; - - is determined  from the IDL "size" function.
;
sdfstr="SDF format"
;hdrinit=2000LL
safety=100LL
np=N_params()
if(np ne 3) then begin
   Message, 'Usage:',/info
   Message, 'sdf_write,fname,label,data,/quiet,hinitsz=hinitsz', /info
   Message, 'Purpose: Write out 1 dataset to an SDF file ', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         label = brief string descriptor of dataset', /info
   Message, '         data  = the data to be written', /info
   Message, '         /quiet  = keyword to suppress non-error messages', /info
   Message, '         hinitsz  -- keyword set to initial header size', /info
   Message, '         (if not set, uses default value of 2000 bytes)', /info
   return
endif

if(keyword_set(quiet)) then begin
   verbose = 0
endif else begin
   verbose = 1
endelse

if(keyword_set(hinitsz)) then begin
   hdrinit=long64(hinitsz)
endif else begin
   hdrinit=2000LL
endelse

get_lun,unit
openr,unit,fname,/swap_if_little_endian,error=error_diag
new_file=0
if(error_diag ne 0) then begin
   if(verbose eq 1) then begin
      print, "sdf_write: file probably doesn't exist, initializing"
   endif
   new_file=1
   hdrsize=long64(hdrinit)
   close,unit
endif

if (new_file eq 0) then begin
   sdftestbyte=bytarr(11)
   readu,unit,sdftestbyte
   sdftest=string(sdftestbyte)

   if( new_file eq 0) then begin
      if(not strcmp(sdftest,sdfstr,10)) then begin
         print,"sdf_write: fname = ",fname," is not an SDF file, returning"
         close,unit
         free_lun,unit
         return
      endif
   endif
   close,unit
endif

if (new_file eq 1) then begin ; write out initial header data positions, etc
   openw,unit,fname,/swap_if_little_endian
   init_bytes=bytarr(11)
   init_bytes[0:9]=byte(sdfstr)
   init_bytes[10]=byte(0)
   point_lun,unit,long64(0)
   writeu,unit,init_bytes
   point_lun,-unit,curpos
   curpos=long64(curpos)
   hdrpos=long64(curpos)+long64(2)*long64(8)+long64(4)+long64(8)
   writeu,unit,hdrpos
   filler=bytarr(hdrsize)
   writeu,unit,filler
   point_lun,-unit,curpos
   datapos=long64(curpos)
   point_lun,unit,hdrpos-long64(8)-long(4)-long64(8)
   writeu,unit,datapos
   next_order=0L
   writeu,unit,next_order
   hdrsize=hdrinit
   writeu,unit,hdrsize
   close,unit
endif

hdrpos=0LL
datapos=0LL
order=0L
hdrsize=0LL

; - get current values of header position, data postion, and dataset order:

openu,unit,fname,/swap_if_little_endian
point_lun,unit,long64(11)
readu,unit,hdrpos
readu,unit,datapos ; datapos is also the current filesize
readu,unit,order
readu,unit,hdrsize
;
; - - check to see if header almost full, if so, print warning
;
if(hdrpos ge (hdrsize-safety)) then begin
   if(verbose eq 1) then begin 
      print,"sdf_write:  header almost full, incrementing hdrsize by hdrinit"
   endif
   oldhdrsize=hdrsize
   hdrsize=hdrsize+hdrinit
   point_lun,unit,long64(strlen(sdfstr)+1+2*8+4)
   writeu,unit,hdrsize
   startloc=long64(strlen(sdfstr)+1+8+oldhdrsize)
   move_up,unit,hdrinit,startloc,datapos
   datapos=datapos+hdrinit
endif
;
;
; - - Now, get dimension information on data:
;
dsize=size(data)
ndim=dsize[0]
if(ndim eq 0) then begin
   dims=lon64arr(1)
   dims[0]=1
   ndim=1
endif else begin
   dims=lon64arr(ndim)
   for ii=0L,ndim-1L do begin
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
endif else if (typecode eq 9) then begin
   dtype = "c"
   nbpw = 8
endif else if ((typecode eq 14) or (typecode eq 15)) then begin
   dtype = "i"
   nbpw = 8
endif else if ((typecode eq 2) or (typecode eq 12)) then begin
   dtype = "i"
   nbpw = 2
endif else if(typecode eq 7) then begin
   print,"string data not allowed in sdf_write.  Convert to byte data"
   close,unit
   free_lun,unit
   return
endif else begin
   print, "illegal data type in sdf_write, typecode = ",typecode
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
   print,"dtype seems blank, returning"
   close,unit
   free_lun,unit
   return
endif
;
; - - make sure that label is blank-free:
;
label=strcompress(label,/remove_all)
;
; - - make sure label has *something* in it:
;
if(strlen(label) eq 0) then begin
   label="NA"
endif
;
; - - create the output string for header
;
data_to_string,outputs,order,label,dtype,nbpw,ndim,dims
;
; - - convert to byte array:
;
hdrlen=strlen(outputs)
hdrbytearr=bytarr(hdrlen+1)
hdrbytearr[0:hdrlen-1]=byte(outputs)
hdrbytearr[hdrlen]=10b ; add what I hope is a line-feed character at end
;
; - - write out the header string
;
point_lun,unit,hdrpos
writeu,unit,hdrbytearr
point_lun,-unit,hdrpos ; get new value of hdrpos
hdrpos=long64(hdrpos)
;
; - - Go to datapos and write data:
;
point_lun,unit,datapos
writeu,unit,data
;
; - -  Update datapos:
;
point_lun,-unit,datapos
datapos=long64(datapos)
;
; - - update order:
;
order=long(order)+1L
;
; - - Update new values of hdrpos,datapos, and order in the file:
;
point_lun,unit,long64(strlen(sdfstr)+1)
writeu,unit,long64(hdrpos)
writeu,unit,long64(datapos)
writeu,unit,long(order)
writeu,unit,long64(hdrsize)
;
; - - we're done; close file and return
;
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
