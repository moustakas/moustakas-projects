pro hff_to_maggies, cat, maggies, ivar, filterlist=filterlist, $
  nJy=nJy, usehawki=usehawki, usemoircs=usemoircs
; jm13nov08siena 

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'hff_to_maggies'
       return
    endif

    filterlist = hff_filterlist(short_filter=filt,/useirac,usehawki=usehawki,$
      usemoircs=usemoircs)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    errtags = filt+'_fluxerr'

; conversion factor from nJy or uJy to AB maggies;
; print, (23.0+9.0)/0.4-48.6
;     31.4000
    if keyword_set(nJy) then fact = 10D^(-0.4D*31.4) else $
      fact = 10D^(-0.4D*23.9)

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       if ftag ne -1 and utag ne -1 then begin
          good = where(cat.(utag) gt 0.0 and cat.(utag) lt 1E3,ngood)
;         print, filt[ib], cat[good[0]].(ftag); *fact
          if (ngood ne 0L) then begin
             maggies[ib,good] = cat[good].(ftag)*fact
             ivar[ib,good] = 1D/(cat[good].(utag)*fact)^2.0
          endif
       endif
    endfor
    
return   
end
