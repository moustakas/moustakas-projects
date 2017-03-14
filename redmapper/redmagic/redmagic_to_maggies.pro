pro redmagic_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
; jm17feb19siena - convert a REDMAGIC-style catalog to maggies

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'redmagic_to_maggies'
       return
    endif

    filterlist = redmagic_filterlist(short=short)
    nbands = n_elements(filterlist)

    tags = 'MODELMAG_'+short
    errtags = 'MODELMAGERR_'+short

; construct maggies and ivarmaggies in each band       
    maggies = fltarr(nbands,ngal)
    ivarmaggies = fltarr(nbands,ngal)
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
