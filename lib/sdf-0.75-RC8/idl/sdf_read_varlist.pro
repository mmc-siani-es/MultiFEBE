; Simple Data Format
; http://solarmuri.ssl.berkeley.edu/~fisher/public/software/SDF
; Copyright (C) 2006,2007 University of California
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
;
PRO sdf_read_varlist,fname,VARS=vars,SUFFIX=suffix

;
; Procedure to read sdf data files 
;   - sdf files use meta data to store variable names/types 
;     so this procedure should work with any sdf file containing datasets
;     with unique (non-repeated) labels
;
ON_ERROR, 2
IF(N_PARAMS() ne 1) THEN BEGIN
   MESSAGE, 'Usage:',/info
   MESSAGE, 'sdf_read_varlist,fname,vars=vars,suffix=suffix', /info
   MESSAGE, 'Purpose: Read data from sdf data file', /info
   MESSAGE, 'Input:   fname  = name of sdf file to read', /info
   MESSAGE, 'Keyword: vars  = string with names of variables to read', /info
   MESSAGE, 'Keyword: suffix  = string to append to variable names', /info
   MESSAGE, 'Output: variables in main IDL session read from file', /info
   MESSAGE, '                                                        ', /info
   MESSAGE, 'Example:  Read all variables in file junk.sdf:',/info
   MESSAGE, 'sdf_read_varlist,"junk.sdf"',/info
   MESSAGE, '                                                       ', /info
   MESSAGE, 'Example:  Read vx and bz from file junk.sdf:',/info
   MESSAGE, 'sdf_read_varlist,"junk.sdf",vars="vx bz"',/info
   RETURN
ENDIF

IF(KEYWORD_SET(suffix)) THEN BEGIN
   sfx = suffix ; keyword if set
ENDIF ELSE BEGIN
   sfx = "" ; null string if not set
ENDELSE

SDF_QUERY, fname, ndat, /quiet ; Get number of datasets in file
IF(KEYWORD_SET(vars)) THEN BEGIN ; find user selected variables
  vartokens = STRSPLIT(vars, " ", /EXTRACT, COUNT=ntoks)
  FOR i=0,ntoks-1 DO BEGIN
     match = SDF_LABMATCH(fname, ndat, vartokens[i])
     order = match[0]
     IF (order NE -1) THEN BEGIN
         SDF_READ, fname, order, label, data, /quiet
         IF(N_ELEMENTS(data) EQ 1) THEN BEGIN
            data = data[0]
         ENDIF
         (SCOPE_VARFETCH(label+sfx, /ENTER, LEVEL=-1)) = data
     ENDIF ELSE BEGIN
       PRINT, vartokens[i] + ' is not a valid variable name'
     ENDELSE
  ENDFOR
ENDIF ELSE BEGIN ; find all variables
  FOR i=0,ndat-1 DO BEGIN
     SDF_READ, fname, i, label, data, /quiet
     IF(N_ELEMENTS(data) EQ 1) THEN BEGIN
        data = data[0]
     ENDIF
     (SCOPE_VARFETCH(label+sfx, /ENTER, LEVEL=-1)) = data
  ENDFOR
ENDELSE
RETURN
END
