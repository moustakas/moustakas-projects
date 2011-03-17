;+
; NAME:
;   BOOTES_TO_MAGGIES
;
; PURPOSE:
;   Convert the BOOTES matched-aperture photometry to maggies. 
;
; INPUTS: 
;   bootes - compatible BOOTES catalog
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   psf - use the PSF magnitudes (default is to use the aperture
;     magnitudes) 
;   nozpoffset - do not apply the derived zeropoint offsets or IRAC
;     aperture corrections
;   itot - use I_tot rather than I_auto (see BUILD_AGES_PHOTOMETRY) 
;
; OUTPUTS: 
;   maggies - 
;   ivarmaggies - 
;
; OPTIONAL OUTPUTS:
;   filterlist - 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jan 05, UCSD
;
; Copyright (C) 2010, John Moustakas
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

pro bootes_to_maggies, bootes, maggies, ivar, psf=psf, use_aper=use_aper, $
  totalmag=totalmag, itot=itot, filterlist=filterlist

    ngal = n_elements(bootes)    
    if (ngal le 0L) then begin
       doc_library, 'bootes_to_maggies'
       return
    endif

    filterlist = bootes_filterlist()
    vega2ab = bootes_vega2ab()
;   zpoffset = bootes_zpoffset(nozpoffset=nozpoffset)
    nbands = n_elements(filterlist)

; correct for Galactic extinction    
    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell,/silent)
    glactc, bootes.i_alpha_j2000, bootes.i_delta_j2000, $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)

; if /TOTALMAG then use the *measured* I-band mag_auto magnitude,
; unless /ITOT is set, in which case use the "kludge" I-band magnitude
; based on the combination of I_Kron and I_R,Kron (see
; BUILD_AGES_PHOTOMETRY); if TOTALMAG is *not* set, then check for the
; PSF keyword, otherwise use the aperture magnitudes (default: 4"
; diameter aperture magnitude)
    bands = ['U','Bw','R','I','z','J','H','Ks','ch1','ch2','ch3','ch4']
    if (n_elements(use_aper) eq 0) then use_aper = '04'
    if keyword_set(totalmag) then begin
       splog, 'Computing total magnitudes via '+use_aper+$
         ' arcsec diameter aperture magnitudes'
       if keyword_set(itot) then itagname = 'i_tot' else itagname = 'i_mag_auto'
       itottag = tag_indx(bootes[0],itagname)
       if (itottag eq -1) then message, 'I-band magnitude '+itagname+' not found!'
       iapertag = tag_indx(bootes[0],'I_mag_aper_'+use_aper)
       iapererrtag = tag_indx(bootes[0],'I_magerr_aper_'+use_aper)
       tags = bands+'_mag_aper_'+use_aper
       errtags = bands+'_magerr_aper_'+use_aper
    endif else begin
       if keyword_set(psf) then begin
          splog, 'Computing PSF magnitudes'
          tags = bands+'_mag_psf'
          errtags = bands+'_magerr_psf'
       endif else begin
          splog, 'Computing '+use_aper+' arcsec diameter aperture magnitudes'
          tags = bands+'_mag_aper_'+use_aper
          errtags = bands+'_magerr_aper_'+use_aper
       endelse
    endelse

; construct maggies and ivarmaggies in each band       
    maggies = fltarr(nbands,ngal)
    ivar = fltarr(nbands,ngal)
    for ii = 0, nbands-1 do begin
       ftag = tag_indx(bootes[0],tags[ii])
       utag = tag_indx(bootes[0],errtags[ii])
       if (ftag[0] eq -1) or (utag[0] eq -1) then begin
          splog, 'No photometry found for band '+tags[ii]+'!'
       endif else begin
          if keyword_set(totalmag) then begin ; total magnitude
             good = where($
               (bootes.(itottag) gt 0.0) and (bootes.(itottag) lt 90.0) and $
               (bootes.(iapertag) gt 0.0) and (bootes.(iapertag) lt 90.0) and $
               (bootes.(iapererrtag) gt 0.0) and (bootes.(iapererrtag) lt 90.0) and $
               (bootes.(ftag) gt 0.0) and (bootes.(ftag) lt 90.0) and $
               (bootes.(utag) gt 0.0) and (bootes.(utag) lt 90.0),ngood)
          endif else begin      ; PSF or aperture magnitude
             good = where($
               (bootes.(ftag) gt 0.0) and (bootes.(ftag) lt 90.0) and $
               (bootes.(utag) gt 0.0) and (bootes.(utag) lt 90.0),ngood)
          endelse
          if (ngood ne 0L) then begin
             magerr = bootes[good].(utag)
             if keyword_set(totalmag) then begin
                mag = bootes[good].(itottag) + (bootes[good].(ftag)-bootes[good].(iapertag)) + $
                  vega2ab[ii] - kl[ii]*ebv[good]; + zpoffset[ii]
             endif else begin
                mag = bootes[good].(ftag) + vega2ab[ii] - kl[ii]*ebv[good]; + zpoffset[ii] 
             endelse
             maggies[ii,good] = 10.0^(-0.4*mag)
             notzero = where((maggies[ii,good] gt 0.0),nnotzero)
             if (nnotzero ne 0L) then ivar[ii,good[notzero]] = $
               1.0/(0.4*alog(10.0)*(maggies[ii,good[notzero]]*magerr[notzero]))^2
          endif
       endelse
    endfor

; apply a minimum photometric error
    minerr = bootes_minerror()
    if (n_elements(minerr) ne nbands) then message, 'Problem here!'
    k_minerror, maggies, ivar, minerr

return   
end
