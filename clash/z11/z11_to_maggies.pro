pro z11_to_maggies, cat, maggies, ivar, filterlist=filterlist
; jm12aug14siena

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'cat_to_maggies'
       return
    endif

    filterlist = z11_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    errtags = filt+'_fluxerr'

; conversion factor from nJy to AB mags
    fact = 10D^(-0.4D*23.9)

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       good = where(cat.(utag) gt 0.0,ngood)
       if (ngood ne 0L) then begin
          maggies[ib,good] = cat[good].(ftag)*fact
          ivar[ib,good] = 1D/(cat[good].(utag)*fact)^2.0
       endif
    endfor
    
return   
end
