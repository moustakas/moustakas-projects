pro alpha_make_softlinks, basepath
; jm08jul07nyu - generate soft links to the raw and pre-processed
; data; only deal about the red side for now; called by
; ALPHA_REDUCE_ALL 

    if (n_elements(basepath) eq 0L) then begin
       splog, 'Please specify BASEPATH'
       return
    endif
    
    if (file_test(basepath+'preproc',/dir) eq 0L) or $
      (file_test(basepath+'rawdata',/dir) eq 0L) or $
      (file_test(basepath+'Raw',/dir) eq 0L) then $
        message, 'Go forth and make the requisite directories'

; assumes [r,b] prefix for RAWLIST and [p] prefix for PROCLIST

    pushd, basepath
    rawlist = file_search('rawdata/r????.fits.gz',count=nraw)    
    proclist = file_search('preproc/p????.fits.gz',count=nproc)
    if (nraw eq 0L) and (nproc eq 0L) then begin
       splog, 'No files found!'
       stop
    endif

    finallist = rawlist
    if (nproc gt 0L) then begin
       these = where_array(strmid(file_basename(proclist),1),$
         strmid(file_basename(rawlist),1))    
       finallist[these] = proclist
    endif

    outlist = 'Raw/'+repstr(file_basename(finallist),'p','r')
;   niceprint, finallist, outlist

    for ii = 0L, nraw-1L do begin
       splog, 'Building '+finallist[ii]+' --> '+outlist[ii]
       spawn, 'ln -sf ../'+finallist[ii]+' '+outlist[ii]
    endfor
    popd
    
return
end
    
