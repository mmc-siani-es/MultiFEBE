pro sdf_expand_size,data,match,data_ex
;
; - - Purpose of this procedure is to create a new array, that is
; - - based on an existing
; - - array, data, and an array of matching datasets, match, with one
; - - additional dimension, corresponding to the number of non-trivial values
; - - of the match array.  The match array is created with the sdf_labmatch
; - - function.
;
; - - First, figure out the parameters (size and type) of the dataset data:

dsize=long64(size(data))
dslen=n_elements(dsize)
dsize_ext=lon64arr(dslen+1)
;
; - - Now, figure out how many non-trivial elements of match there are:
;
nontrivial=where(match ne -1)
nsig=n_elements(nontrivial)
if(nsig le 0L) then begin
   print,"sdf_expand_size: error, no non-trivial matches, no expanded array"
   return
endif
nontrivial=nontrivial[0:nsig-1]
index=match[nontrivial]
ndimdat=dsize[0]
dsize_ext[0]=ndimdat+1LL
for i=1L,ndimdat do begin
  dsize_ext[i]=dsize[i]
endfor
dsize_ext[ndimdat+1]=long64(nsig)
for i=ndimdat+1,dslen-1 do begin
  dsize_ext[i+1]=dsize[i]
endfor
dsize_ext[dslen]=dsize_ext[dslen]*nsig
data_ex=make_array(size=dsize_ext);
end
