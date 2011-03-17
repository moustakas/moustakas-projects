;+
; NAME:
;   MZ_TO_MAGGIES
;
; PURPOSE:
;   Convert the BOOTES matched-aperture photometry to maggies
;   specifically for the AGES/MZ project.
;
; INPUTS: 
;   bootes - BOOTES photometric catalog
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
;   J. Moustakas, 2010 Oct 18, UCSD
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

pro mz_to_maggies, cat, maggies, ivar, psf=psf, use_aper=use_aper, $
  totalmag=totalmag, itot=itot, filterlist=filterlist

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'mz_to_maggies'
       return
    endif

    filterlist = mz_filterlist(bands=bands,vega2ab=vega2ab,$
      zpoffset=zpoffset,minerr=minerr,nozpoffset=nozpoffset)
    nfilt = n_elements(filterlist)

; correct for Galactic extinction    
    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell,/silent)
    glactc, cat.i_alpha_j2000, cat.i_delta_j2000, $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)

; if /TOTALMAG then use the *measured* I-band mag_auto magnitude,
; unless /ITOT is set, in which case use the "kludge" I-band magnitude
; based on the combination of I_Kron and I_R,Kron (see
; BUILD_AGES_PHOTOMETRY); if TOTALMAG is *not* set, then check for the
; PSF keyword, otherwise use the aperture magnitudes (default: 4"
; diameter aperture magnitude)
    if (n_elements(use_aper) eq 0) then use_aper = '04'
    if keyword_set(totalmag) then begin
       splog, 'Computing total magnitudes via '+use_aper+$
         ' arcsec diameter aperture magnitudes'
       if keyword_set(itot) then itagname = 'i_tot' else itagname = 'i_mag_auto'
       itottag = tag_indx(cat[0],itagname)
       if (itottag eq -1) then message, 'I-band magnitude '+itagname+' not found!'
       iapertag = tag_indx(cat[0],'I_mag_aper_'+use_aper)
       iapererrtag = tag_indx(cat[0],'I_magerr_aper_'+use_aper)
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
    maggies = fltarr(nfilt,ngal)
    ivar = fltarr(nfilt,ngal)
    for ii = 0, nfilt-1 do begin
       ftag = tag_indx(cat[0],tags[ii])
       utag = tag_indx(cat[0],errtags[ii])
       if (ftag[0] eq -1) or (utag[0] eq -1) then begin
          splog, 'No photometry found for band '+tags[ii]+'!'
       endif else begin
          if keyword_set(totalmag) then begin ; total magnitude
             good = where($
               (cat.(itottag) gt 0.0) and (cat.(itottag) lt 90.0) and $
               (cat.(iapertag) gt 0.0) and (cat.(iapertag) lt 90.0) and $
               (cat.(iapererrtag) gt 0.0) and (cat.(iapererrtag) lt 90.0) and $
               (cat.(ftag) gt 0.0) and (cat.(ftag) lt 90.0) and $
               (cat.(utag) gt 0.0) and (cat.(utag) lt 90.0),ngood)
          endif else begin      ; PSF or aperture magnitude
             good = where($
               (cat.(ftag) gt 0.0) and (cat.(ftag) lt 90.0) and $
               (cat.(utag) gt 0.0) and (cat.(utag) lt 90.0),ngood)
          endelse
          if (ngood ne 0L) then begin
             magerr = cat[good].(utag)
             if keyword_set(totalmag) then begin
                mag = cat[good].(itottag) + (cat[good].(ftag)-cat[good].(iapertag)) + $
                  vega2ab[ii] - kl[ii]*ebv[good] + zpoffset[ii]
             endif else begin
                mag = cat[good].(ftag) + vega2ab[ii] - kl[ii]*ebv[good] + zpoffset[ii] 
             endelse
             maggies[ii,good] = 10.0^(-0.4*mag)
             notzero = where((maggies[ii,good] gt 0.0),nnotzero)
             if (nnotzero ne 0L) then ivar[ii,good[notzero]] = $
               1.0/(0.4*alog(10.0)*(maggies[ii,good[notzero]]*magerr[notzero]))^2
          endif
       endelse
    endfor

; the z-band photometry of objects near the edges of the individual
; mosaics are suspect, so toss them out here
    zfile = ages_path(/mycat)+'zbootes/zbootes_field_centers.fits.gz'
    splog, 'Reading '+zfile
    field = mrdfits(zfile,1)
    nfield = n_elements(field)
    zmask = intarr(ngal) ; all start out bad
    pad = 10.0/3600.0 ; 10" padding
    
    for ii = 0, nfield-1 do zmask = zmask or $
      ((cat.ra gt field[ii].ra-field[ii].dra/2.0+pad) and $
      (cat.ra lt field[ii].ra+field[ii].dra/2.0-pad) and $
      (cat.dec gt field[ii].dec-field[ii].ddec/2.0+pad) and $      
      (cat.dec lt field[ii].dec+field[ii].ddec/2.0-pad))
;   good = where(zmask);,comp=bad)
    zfilt = (where(strmatch(filterlist,'*bok*')))[0]
    bad = where((zmask eq 0) and (maggies[zfilt,*] gt 0.0))

;   djs_plot, cat.ra, cat.dec, ps=3, ysty=3               
;   djs_oplot, cat[good].ra, cat[good].dec, ps=3, color='blue'
;   djs_oplot, cat[bad].ra, cat[bad].dec, ps=3, color='red'
    maggies[zfilt,bad] = 0.0
    ivar[zfilt,bad] = 0.0
    
; apply a minimum photometric error
    k_minerror, maggies, ivar, minerr

return   
end
