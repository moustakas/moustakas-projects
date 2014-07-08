pro alhambra_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
; jm14jul07siena - convert an ALHAMBRA-style catalog to maggies

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'alhambra_to_maggies'
       return
    endif

    filterlist = alhambra_filterlist(short=short)
    nbands = n_elements(filterlist)

    tags = short
    errtags = 'd'+short

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivarmaggies = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       good = where(cat.(ftag) gt 0.0 and cat.(ftag) lt 90 and $
         cat.(utag) gt 0.0 and cat.(utag) lt 90,ngood)
       if (ngood ne 0L) then begin
          maggies[ib,good] = 10D^(-0.4D*cat[good].(ftag))
          magerr = cat[good].(utag)
          ivarmaggies[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
       endif
    endfor
    
return   
end
