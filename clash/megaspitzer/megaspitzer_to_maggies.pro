pro megaspitzer_to_maggies, clash, maggies, ivar, filterlist=filterlist, $
  nominerror=nominerror, zpt=zpt
; jm12jul25siena
    
    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'megaspitzer_to_maggies'
       return
    endif

    filterlist = megaspitzer_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    errtags = filt+'_fluxerr'

; conversion factor from nJy to AB mags
    fact = 10D^(-0.4D*31.4)

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(clash[0],tags[ib])
       utag = tag_indx(clash[0],errtags[ib])

       good = where(clash.(utag) gt 0.0,ngood)
       if (ngood ne 0L) then begin
          maggies[ib,good] = clash[good].(ftag)*fact
          ivar[ib,good] = 1D/(clash[good].(utag)*fact)^2.0
       endif
    endfor
    
; apply a minimum photometric error
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.1,nbands)
       isirac = where(strmatch(filterlist,'*irac*',/fold))
       if isirac[0] ne -1 then minerr[isirac] = 0.0 ; note!
       k_minerror, maggies, ivar, minerr
    endif
    
return   
end
