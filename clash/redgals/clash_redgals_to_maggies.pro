pro clash_redgals_to_maggies, cat, maggies, ivar, filterlist=filterlist
; jm15sep02siena - fluxes have been corrected for dust attenuation

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'clash_redgals_to_maggies'
       return
    endif

    filterlist = clash_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    ivartags = filt+'_ivar'

; conversion factor from uJy to AB maggies
    fact = 10D^(-0.4D*23.9)

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],ivartags[ib])

       if ftag ne -1 and utag ne -1 then begin
          good = where(cat.(utag) gt 0.0,ngood)
;         good = where(cat.(utag) gt 0.0 and cat.(utag) lt 1E3,ngood)
          if (ngood ne 0L) then begin
             maggies[ib,good] = float(cat[good].(ftag)*fact)
             ivar[ib,good] = float(cat[good].(utag)/fact^2.0)
          endif
       endif
    endfor

;; minimum photometric error    
;    minerr = fltarr(nbands)+0.02
;    k_minerror, maggies, ivar, minerr
    
return   
end
