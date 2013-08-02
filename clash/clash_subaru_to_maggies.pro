pro clash_subaru_to_maggies, cat, maggies, ivar, filterlist=filterlist, $
  nominerror=nominerror
; jm13jul17siena - convert the CLASH/Subaru photometry to maggies 

    ngal = n_elements(cat)
    if ngal eq 0L then return

    filterlist = clash_subaru_filterlist()
    nbands = n_elements(filterlist)

    tags = ['B','V','R','I','Z']
    errtags = 'd'+['B','V','R','I','Z']

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       good = where(cat.(utag) gt 0.0,ngood)
       if (ngood ne 0L) then begin
          magerr = cat[good].(utag)
          maggies[ib,good] = 10.0^(-0.4*cat[good].(ftag))
          ivar[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
       endif
    endfor
    
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.02,nbands)
       k_minerror, maggies, ivar, minerr
    endif

return   
end


