pro ndwfs_to_maggies, ndwfs, maggies, ivar, filterlist=filterlist, $
  nozpoffset=nozpoffset
; jm09aug25ucsd - written

    ngal = n_elements(ndwfs)    
    if (ngal le 0L) then begin
       doc_library, 'ndwfs_to_maggies'
       return
    endif

    filterlist = ndwfs_filterlist()
    nbands = n_elements(filterlist)

    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell)
    glactc, ndwfs.i_alpha_j2000, ndwfs.i_delta_j2000, $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)

    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
    zpoffset = ndwfs_zpoffset()
    if keyword_set(nozpoffset) then zpoffset = zpoffset*0.0
;   print, zpoffset
    
    maggies = fltarr(nbands,ngal)
    ivar = fltarr(nbands,ngal)

    tags = ['Bw','R','I','K']+'_MAG_AUTO'
    errtags = ['Bw','R','I','K']+'_MAGERR_AUTO'

    for ii = 0L, nbands-1L do begin
       ftag = tag_indx(ndwfs[0],tags[ii])
       utag = tag_indx(ndwfs[0],errtags[ii])
       if (ftag[0] eq -1) or (utag[0] eq -1) then begin
          splog, 'No photometry found for band '+tags[ii]+'!'
       endif else begin
          good = where($
            (ndwfs.(ftag) gt 0.0) and (ndwfs.(ftag) lt 90.0) and $
            (ndwfs.(utag) gt 0.0) and (ndwfs.(utag) lt 90.0),ngood)
          if (ngood ne 0L) then begin
             mag = ndwfs[good].(ftag) - kl[ii]*ebv[good] + vega2ab[ii] + zpoffset[ii]
             magerr = ndwfs[good].(utag)
             maggies[ii,good] = 10.0^(-0.4*mag)
             notzero = where((maggies[ii,good] gt 0.0),nnotzero)
             if (nnotzero ne 0L) then ivar[ii,good[notzero]] = $
               1.0/(0.4*alog(10.0)*(maggies[ii,good[notzero]]*magerr[notzero]))^2
          endif
       endelse
    endfor

; apply a minimum photometric error
    minerr = [0.05,0.05,0.05,0.10]
    primus_minerror, maggies, ivar, minerr

return   
end
