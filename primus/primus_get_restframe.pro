function primus_get_restframe, slits
; jm08jun15nyu - given a PRIMUS slits structure, this returns the
;                rest-frame (k-corrected) data structure (currently
;                only works for DEEP2)

    nslits = n_elements(slits)
    if (nslits eq 0L) then begin
       doc_library, 'primus_get_restframe'
       return, -1L
    endif
    
    alldeep = read_deep2(/kcorr)
    rest = im_empty_structure(alldeep[0],ncopies=nslits,$
      empty_value=-999.0)

; DEEP2
    
    these = where(slits.deep2_zcat,nthese)
    if (nthese eq 0L) then return, rest

; cross-match the PRIMUS and DEEP2 catalogs; there may be
; repeats, so don't use SPHEREMATCH here    
    
    pivot = min([slits[these].objno,alldeep.objno])
    maxval = max([slits[these].objno,alldeep.objno]-pivot)+1L
    thismatch = lonarr(maxval)-1L
    thismatch[alldeep.objno-pivot] = lindgen(n_elements(alldeep))
    match = thismatch[slits[these].objno-pivot]
    
; slow-ass FOR loop; MATCH2 is not very intuitive
;   match = lonarr(nthese)
;   for ii = 0L, nthese-1L do match[ii] = $
;     where(slits[these[ii]].objno eq alldeep.objno)
;   match2, slits[these].objno, alldeep.objno, m1, m2

    flag = where(match eq -1L,nflag,comp=good,ncomp=ngood)
    if (nflag ne 0L) then begin
       for iflag = 0L, n_elements(flag)-1L do $
         splog, 'Object(s) '+string(slits[these[flag[iflag]]].objno,$
         format='(I0)')+' not found in my DEEP2 catalog!'
    endif
    if (ngood eq 0L) then message, 'This should not happen!' else begin
       these = these[good]
       match = match[good]
    endelse

    rest[these] = alldeep[match]

return, rest
end
