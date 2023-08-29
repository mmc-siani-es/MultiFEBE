pro move_dn,unit,shiftdn,startloc,endloc
; move the contents of a file down shiftdn bytes starting at byte startloc
; shiftdn, startloc, endloc should be 64 bit integers
bigbuf=1048576LL
sz_shift=long64(endloc-startloc-shiftdn)
nbuf=sz_shift / bigbuf
fag_end = sz_shift mod bigbuf
lower_loc=startloc
point_lun,unit,lower_loc
if(sz_shift gt 0) then begin
  for i=0LL,nbuf do begin
     bufsize = (fag_end > ((1LL < (nbuf-i))*bigbuf))
     if(bufsize gt 0) then begin
        buf=bytarr(bufsize)
;       print,"bufsize = ",bufsize
        point_lun,unit,lower_loc+shiftdn
        readu,unit,buf
        point_lun,unit,lower_loc
        writeu,unit,buf
        lower_loc=lower_loc+bufsize
     endif
  endfor
  point_lun,unit,endloc-shiftdn
endif
truncate_lun,unit
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
