pro zbootes_to_maggies, zbootes, maggies, ivar, filterlist=filterlist, $
  nozpoffset=nozpoffset
; jm09aug25ucsd - written

    ngal = n_elements(zbootes)    
    if (ngal le 0L) then begin
       doc_library, 'zbootes_to_maggies'
       return
    endif

    filterlist = zbootes_filterlist()
    nbands = n_elements(filterlist)

    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell)
    glactc, zbootes.alpha_j2000, zbootes.delta_j2000, $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)

    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)*0.0 ; already in AB
    zpoffset = 0.0
    if keyword_set(nozpoffset) then zpoffset = zpoffset*0.0
    
    maggies = fltarr(nbands,ngal)
    ivar = fltarr(nbands,ngal)

    tags = 'MAG_AUTO'
    errtags = 'MAGERR_AUTO'

; see Cool+07 for the quality cuts used    
    for ii = 0L, nbands-1L do begin
       ftag = tag_indx(zbootes[0],tags[ii])
       utag = tag_indx(zbootes[0],errtags[ii])
       good = where((zbootes.nobs gt 5) and (zbootes.photflag eq 0) and $
         (zbootes.(ftag) gt 0.0) and (zbootes.(ftag) lt 90.0) and $
         (zbootes.(utag) gt 0.0) and (zbootes.(utag) lt 90.0),ngood)
       if (ngood ne 0L) then begin
          mag = zbootes[good].(ftag) - kl[ii]*ebv[good] + vega2ab[ii] + zpoffset[ii]
          magerr = zbootes[good].(utag)
          maggies[ii,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[ii,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivar[ii,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[ii,good[notzero]]*magerr[notzero]))^2
       endif
    endfor

; apply a minimum photometric error
    minerr = 0.05
    primus_minerror, maggies, ivar, minerr

return   
end
