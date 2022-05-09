function sdf_read_arr,fname,labmatch
;
; - - Reads in all datasets with a common value of label, and creates
; - - a new array, with the number of dimensions of data_arr 1 greater than the
; - - dimensions of the original data written out. The array data_arr
; - - correponds to all occurrences of data for which label matches labmatch.
; - - This is intended to create time-series of quantities which
; - - are written repeatedly in an sdf file.  For example, in a thin flux
; - - tube simulation, where the scalar loop length is written out at each
; - - output time step, using sdf_read_arr will result in an array of all
; - - values of loop length over time.  If data is multiply dimensioned to begin
; - - with, one additional dimension is added as the last dimension.  The
; - - value of this last dimension is the number of times the selected dataset
; - - label appears in the sdf file.
;

np=N_params()
if(np ne 2) then begin
   Message, 'Usage:',/info
   Message, 'data=sdf_read_arr(fname,labtest)', /info
   Message, 'Purpose: Read *all* datasets from file for which labtest', /info
   Message, '         matches dataset label', /info
   Message, 'Input:   fname  = the name of the SDF file', /info
   Message, '         labtest = input label for desired dataset in file', /info
   Message, 'Output:  data_arr = the returned data, as a "time series"  ', /info
   Message, 'Note:    data_arr has 1 more dimension than the individual', /info
   Message, '         datasets. The value of the last dimension is the', /info
   Message, '         number of occurrences of label in the sdf file.', /info
   return, -1
endif

match=sdf_labmatch(fname,ndat,labmatch)

;
; - - Test to make sure that there's at least 1 occurrence of label
;

if(match[0] eq -1L) then begin
  print,"sdf_read_arr:  No label matches, returning -1 to data_arr"
  data_arr=-1L
  return, data_arr
endif

nontrivial=where(match ne -1)
nsig=long(n_elements(nontrivial))
nontrivial=nontrivial[0:nsig-1]
index=long(match[nontrivial])

sdf_read,fname,index[0],label,data,/quiet
if(nsig eq 1L) then begin
   print,"sdf_read_arr:  Only one match.  data_arr is single matching dataset"
   data_arr=data
   return, data_arr
endif

sdf_expand_size,data,index,data_arr

;
; - - To read in the recurring datasets without knowing in advance how
; - - many dimensions they have, use the n_elements function and revert
; - - to 1-d indexing, assuming column major memory mapping for IDL.
;

offset=n_elements(data)

for i=0,nsig-1L do begin
   sdf_read,fname,index[i],label,data,/quiet
   data_arr[i*offset:(i+1)*offset-1]=data[0:offset-1]
endfor

;
; - - in case data_arr is a 1 x X array, just make it an X dim'l 1d array
;

data_arr=reform(data_arr)

return, data_arr
end
