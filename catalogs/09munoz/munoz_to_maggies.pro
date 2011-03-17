;+
; NAME:
;   MUNOZ_TO_MAGGIES
;
; PURPOSE:
;   Convert the Munoz-Mateos et al. (2009) photometry to
;   extinction-corrected AB maggies. 
;
; INPUTS: 
;   cat - input photometric catalog [NGAL] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   maggies - output maggies [10,NGAL]
;   ivarmaggies - corresponding inverse variance array [10,NGAL]  
;   filterlist - filter set
;
; COMMENTS:
;   A minimum error of 0.02 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Oct 30 - written
;
; Copyright (C) 2009, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro munoz_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

    nobj = n_elements(cat)
    if (nobj eq 0) then begin
       doc_library, 'munoz_to_maggies'
       return
    endif
    
    names = ['fuv','nuv','u','g','r','i','z','j','h','ks']
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      twomass_filterlist()]
    minerrors = replicate(0.05,n_elements(names))

; zeropoint and AB corrections; photometry has already been corrected
; for foreground Galactic extinction
    zpt = minerrors*0.0
    vega2ab = minerrors*0.0 ; already in AB

;   weff = k_lambda_eff(filterlist=filterlist)
;   kl = k_lambda(weff,/odonnell,R_V=3.1)
;   glactc, cat.ra, cat.dec, 2000.0, gl, gb, 1, /deg
;   ebv = dust_getval(gl,gb,/interp,/noloop)
    ebv = fltarr(nobj)

    nband = n_elements(names)
    maggies = dblarr(nband,nobj)
    ivarmaggies = dblarr(nband,nobj)
    for iband = 0L, nband-1 do begin
       ftag = tag_indx(cat[0],names[iband]+'')
       utag = tag_indx(cat[0],names[iband]+'_err')
       good = where((cat.(ftag) gt 0.0) and $
         (cat.(utag) gt 0.0) and (cat.(ftag) lt 90.0) and $
         (cat.(utag) lt 90.0),ngood)
       if (ngood ne 0) then begin
          mag = cat[good].(ftag)-kl[iband]*ebv[good]+vega2ab[iband]+zpt[iband]
          magerr = cat[good].(utag)
          maggies[iband,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[iband,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,good[notzero]]*magerr[notzero]))^2
       endif 
    endfor

    k_minerror, maggies, ivarmaggies, minerrors

return
end
