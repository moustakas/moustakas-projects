;+
; NAME:
;   PARSE_NDWFS_PHOT()
;
; PURPOSE:
;   Compute aperture colors and total magnitudes from the NDWFS BwRIK 
;   catalogs.
;
; INPUTS: 
;   bwband, rband, iband, kband - AGES/NDWFS (FITS-format) catalogs
;     (e.g., ages_path(/analysis)+'catalog.ndwfsb.fits.gz')
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   allndwfs - by default this routine takes the four *separate* NDWFS
;     catalogs pertaining to the AGES survey; however, this routine
;     can also be called with /ALLNDWFS in order to compute total
;     magnitudes and colors for the full parent set of NDWFS galaxies
;     (see, specifically, the output from UNPACK_NDWFS, and
;     NDWFS_KBAND_TESTS for an example of how to call this routine);
;     see BUILD_AGES_PHOTOMETRY for the default usage
;
; OUTPUTS: 
;   ndwfs - line-match data structure with lots of goodies
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Sep 14, UCSD - written
;   jm09nov09ucsd - documented
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

function parse_ndwfs_phot, bwband, rband, iband, kband, $
  allndwfs=allndwfs, nozpoffset=nozpoffset

    if keyword_set(allndwfs) then begin
       if (n_params() ne 1) then begin
          splog, 'With /ALLNDWFS pass the output from UNPACK_NDWFS'
          return, -1
       endif
    endif else begin
       if (n_params() ne 4) then begin
          splog, 'BWBAND, RBAND, IBAND, and KBAND input catalogs needed'
          return, -1
       endif
    endelse
    
; NDWFS zeropoint corrections
    zpoffset = ndwfs_zpoffset()
    if keyword_set(nozpoffset) then zpoffset = zpoffset*0.0
    zpbw = zpoffset[0]
    zpr = zpoffset[1]
    zpi = zpoffset[2]
    zpk = zpoffset[3]

    mag_err_floor = 0.001 ; some objects have zero error!!
    
; note: _APER3 = 4", _APER5 = 6" diameter
    
    ngal = n_elements(iband)
    ndwfs = {$
      bw_flags:         -999, $
      r_flags:          -999, $
      i_flags:          -999, $
      k_flags:          -999, $
      
      bw_splitmatch:    -999, $
      r_splitmatch:     -999, $
      i_splitmatch:     -999, $
      k_splitmatch:     -999, $
      
      bwmag_auto:     -999.0, $ ; Kron
      bwmag_auto_err: -999.0, $
      rmag_auto:      -999.0, $
      rmag_auto_err:  -999.0, $
      imag_auto:      -999.0, $
      imag_auto_err:  -999.0, $
      kmag_auto:      -999.0, $
      kmag_auto_err:  -999.0, $

      bwcolor_4:      -999.0, $ ; (Bw-I) color in 4 arcsec aperture
      bwcolor_4_err:  -999.0, $
      rcolor_4:       -999.0, $
      rcolor_4_err:   -999.0, $
      kcolor_4:       -999.0, $
      kcolor_4_err:   -999.0, $

      bwcolor_6:      -999.0, $ ; (Bw-I) color in 6 arcsec aperture
      bwcolor_6_err:  -999.0, $
      rcolor_6:       -999.0, $
      rcolor_6_err:   -999.0, $
      kcolor_6:       -999.0, $
      kcolor_6_err:   -999.0, $

      weight:         -999.0, $
      itot:           -999.0, $ ; total, I-band corrected magnitude
      itot_err:       -999.0, $

      bwtot_4:        -999.0, $ ; total magnitudes using the 4" and 6" aperture colors
      bwtot_4_err:    -999.0, $
      rtot_4:         -999.0, $
      rtot_4_err:     -999.0, $ 
      ktot_4:         -999.0, $
      ktot_4_err:     -999.0, $
      bwtot_6:        -999.0, $
      bwtot_6_err:    -999.0, $
      rtot_6:         -999.0, $
      rtot_6_err:     -999.0, $ 
      ktot_6:         -999.0, $
      ktot_6_err:     -999.0, $

      kcolor_4_r:     -999.0, $ ; for testing
      kcolor_6_r:     -999.0, $
      ktot_r_4:       -999.0, $
      ktot_r_6:       -999.0}

    ndwfs = replicate(ndwfs,ngal)

; SE flags
    good = where(bwband.bw_flags ge 0)
    ndwfs[good].bw_flags = bwband[good].bw_flags

    good = where(rband.r_flags ge 0)
    ndwfs[good].r_flags = rband[good].r_flags

    good = where(iband.i_flags ge 0)
    ndwfs[good].i_flags = iband[good].i_flags

    good = where(kband.k_flags ge 0)
    ndwfs[good].k_flags = kband[good].k_flags

; splitmatch bit
    ndwfs.bw_splitmatch = bwband.bw_flag_splitmatch
    ndwfs.r_splitmatch = rband.r_flag_splitmatch
    ndwfs.i_splitmatch = iband.i_flag_splitmatch
    ndwfs.k_splitmatch = kband.k_flag_splitmatch

; mag_auto
    good = where((bwband.bw_mag_auto gt 0.0) and (bwband.bw_mag_auto lt 90.0) and $
      (bwband.bw_flag_duplicate eq 0)); and (bwband.bw_flag_splitmatch eq 0))
    ndwfs[good].bwmag_auto = bwband[good].bw_mag_auto + zpbw
    ndwfs[good].bwmag_auto_err = bwband[good].bw_magerr_auto > mag_err_floor

    good = where((rband.r_mag_auto gt 0.0) and (rband.r_mag_auto lt 90.0) and $
      (rband.r_flag_duplicate eq 0)); and (rband.r_flag_splitmatch eq 0))
    ndwfs[good].rmag_auto = rband[good].r_mag_auto + zpr
    ndwfs[good].rmag_auto_err = rband[good].r_magerr_auto > mag_err_floor

    good = where((iband.i_mag_auto gt 0.0) and (iband.i_mag_auto lt 90.0) and $
      (iband.i_flag_duplicate eq 0)); and (iband.i_flag_splitmatch eq 0))
    ndwfs[good].imag_auto = iband[good].i_mag_auto + zpi
    ndwfs[good].imag_auto_err = iband[good].i_magerr_auto > mag_err_floor

    good = where((kband.k_mag_auto gt 0.0) and (kband.k_mag_auto lt 90.0) and $
      (kband.k_flag_duplicate eq 0)); and (kband.k_flag_splitmatch eq 0))
    ndwfs[good].kmag_auto = kband[good].k_mag_auto + zpk
    ndwfs[good].kmag_auto_err = kband[good].k_magerr_auto > mag_err_floor

; aperture colors    
    if keyword_set(allndwfs) then begin
       good = where($
         (bwband.bw_flag_duplicate eq 0) and $ ; (bwband.bw_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (bwband.bw_mag_aper_03 gt 0.0) and (bwband.bw_mag_aper_03 lt 90.0) and $
         (bwband.bw_mag_aper_05 gt 0.0) and (bwband.bw_mag_aper_05 lt 90.0) and $
         (iband.i_mag_aper_03 gt 0.0) and (iband.i_mag_aper_03 lt 90.0) and $
         (iband.i_mag_aper_05 gt 0.0) and (iband.i_mag_aper_05 lt 90.0))
       ndwfs[good].bwcolor_4 = bwband[good].bw_mag_aper_03-iband[good].i_mag_aper_03 + (zpbw-zpi)
       ndwfs[good].bwcolor_6 = bwband[good].bw_mag_aper_05-iband[good].i_mag_aper_05 + (zpbw-zpi)
       ndwfs[good].bwcolor_4_err = bwband[good].bw_magerr_aper_03
       ndwfs[good].bwcolor_6_err = bwband[good].bw_magerr_aper_05
    endif else begin
       good = where($
         (bwband.bw_flag_duplicate eq 0) and $ ; (bwband.bw_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (bwband.bw_mag_aper3 gt 0.0) and (bwband.bw_mag_aper3 lt 90.0) and $
         (bwband.bw_mag_aper5 gt 0.0) and (bwband.bw_mag_aper5 lt 90.0) and $
         (iband.i_mag_aper3 gt 0.0) and (iband.i_mag_aper3 lt 90.0) and $
         (iband.i_mag_aper5 gt 0.0) and (iband.i_mag_aper5 lt 90.0))
       ndwfs[good].bwcolor_4 = bwband[good].bw_mag_aper3-iband[good].i_mag_aper3 + (zpbw-zpi)
       ndwfs[good].bwcolor_6 = bwband[good].bw_mag_aper5-iband[good].i_mag_aper5 + (zpbw-zpi)
       ndwfs[good].bwcolor_4_err = bwband[good].bw_magerr_aper3
       ndwfs[good].bwcolor_6_err = bwband[good].bw_magerr_aper5
    endelse

    if keyword_set(allndwfs) then begin
       good = where($
         (rband.r_flag_duplicate eq 0) and $ ; (rband.r_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (rband.r_mag_aper_03 gt 0.0) and (rband.r_mag_aper_03 lt 90.0) and $
         (rband.r_mag_aper_05 gt 0.0) and (rband.r_mag_aper_05 lt 90.0) and $
         (iband.i_mag_aper_03 gt 0.0) and (iband.i_mag_aper_03 lt 90.0) and $
         (iband.i_mag_aper_05 gt 0.0) and (iband.i_mag_aper_05 lt 90.0))
       ndwfs[good].rcolor_4 = rband[good].r_mag_aper_03-iband[good].i_mag_aper_03 + (zpr-zpi)
       ndwfs[good].rcolor_6 = rband[good].r_mag_aper_05-iband[good].i_mag_aper_05 + (zpr-zpi)
       ndwfs[good].rcolor_4_err = rband[good].r_magerr_aper_03
       ndwfs[good].rcolor_6_err = rband[good].r_magerr_aper_05
    endif else begin
       good = where($
         (rband.r_flag_duplicate eq 0) and $ ; (rband.r_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (rband.r_mag_aper3 gt 0.0) and (rband.r_mag_aper3 lt 90.0) and $
         (rband.r_mag_aper5 gt 0.0) and (rband.r_mag_aper5 lt 90.0) and $
         (iband.i_mag_aper3 gt 0.0) and (iband.i_mag_aper3 lt 90.0) and $
         (iband.i_mag_aper5 gt 0.0) and (iband.i_mag_aper5 lt 90.0))
       ndwfs[good].rcolor_4 = rband[good].r_mag_aper3-iband[good].i_mag_aper3 + (zpr-zpi)
       ndwfs[good].rcolor_6 = rband[good].r_mag_aper5-iband[good].i_mag_aper5 + (zpr-zpi)
       ndwfs[good].rcolor_4_err = rband[good].r_magerr_aper3
       ndwfs[good].rcolor_6_err = rband[good].r_magerr_aper5
    endelse       

    if keyword_set(allndwfs) then begin
       good = where($
         (kband.k_flag_duplicate eq 0) and $ ; (kband.k_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (kband.k_mag_aper_03 gt 0.0) and (kband.k_mag_aper_03 lt 90.0) and $
         (kband.k_mag_aper_05 gt 0.0) and (kband.k_mag_aper_05 lt 90.0) and $
         (iband.i_mag_aper_03 gt 0.0) and (iband.i_mag_aper_03 lt 90.0) and $
         (iband.i_mag_aper_05 gt 0.0) and (iband.i_mag_aper_05 lt 90.0))
       ndwfs[good].kcolor_4 = kband[good].k_mag_aper_03-iband[good].i_mag_aper_03 + (zpk-zpi)
       ndwfs[good].kcolor_6 = kband[good].k_mag_aper_05-iband[good].i_mag_aper_05 + (zpk-zpi)
       ndwfs[good].kcolor_4_err = kband[good].k_magerr_aper_03
       ndwfs[good].kcolor_6_err = kband[good].k_magerr_aper_05
    endif else begin
       good = where($
         (kband.k_flag_duplicate eq 0) and $ ; (kband.k_flag_splitmatch eq 0) and $
         (iband.i_flag_duplicate eq 0) and $ ; (iband.i_flag_splitmatch eq 0) and $
         (kband.k_mag_aper3 gt 0.0) and (kband.k_mag_aper3 lt 90.0) and $
         (kband.k_mag_aper5 gt 0.0) and (kband.k_mag_aper5 lt 90.0) and $
         (iband.i_mag_aper3 gt 0.0) and (iband.i_mag_aper3 lt 90.0) and $
         (iband.i_mag_aper5 gt 0.0) and (iband.i_mag_aper5 lt 90.0))
       ndwfs[good].kcolor_4 = kband[good].k_mag_aper3-iband[good].i_mag_aper3 + (zpk-zpi)
       ndwfs[good].kcolor_6 = kband[good].k_mag_aper5-iband[good].i_mag_aper5 + (zpk-zpi)
       ndwfs[good].kcolor_4_err = kband[good].k_magerr_aper3
       ndwfs[good].kcolor_6_err = kband[good].k_magerr_aper5
    endelse

; compute the total I-band magnitude, following Eisenstein+09
;   final = where(ndwfs.imag_auto gt 0.0)
    final = where((ndwfs.rmag_auto gt 0.0) and (ndwfs.rcolor_6 gt -900.0) and $
      (ndwfs.imag_auto gt 0.0))
    
    ikron = ndwfs[final].imag_auto ; I-band Kron
    irband = ndwfs[final].rmag_auto - ndwfs[final].rcolor_6 ; I_R = R_Kron-(R-I)_6"
    iavg = (ikron + irband)/2.0  ; average the two magnitudes
    ifainter = ikron>irband      ; choose the fainter of the two

; compute a statistic that allows us to linearly combine the two
; I-band magnitudes (I_avg and I_fainter); when I_Kron and I_R differ
; by a lot, use I_fainter, otherwise use I_avg    
    weight = exp(-((ikron-irband)/0.2)^2.0)

    itot_flux = (1.0-weight)*10^(-0.4*ifainter) + weight*10^(-0.4*iavg)
    itot = -2.5*alog10(itot_flux)

    ndwfs[final].weight = weight
    ndwfs[final].itot = itot
    ndwfs[final].itot_err = ndwfs[final].imag_auto_err ; conservative error 

; special case for objects with I-band photometry but no R-band
; photometry, in which case the procedure above fails; these are
; presumably in regions where R-band imaging was not obtained
    special = where((ndwfs.rmag_auto lt 0.0) and (ndwfs.imag_auto gt 0.0))
    ndwfs[special].itot = ndwfs[special].imag_auto
    ndwfs[special].itot_err = ndwfs[special].imag_auto_err

; final special case (1 object!) where the R-I color is not defined
; because the 6" R-band aperture color is not defined, but the
; I-band photometry is OK
    special = where((ndwfs.rmag_auto gt 0.0) and (ndwfs.imag_auto gt 0.0) and $
      (ndwfs.rcolor_6 lt 0.0) and (ndwfs.itot lt 0.0))
    ndwfs[special].itot = ndwfs[special].imag_auto
    ndwfs[special].itot_err = ndwfs[special].imag_auto_err
    
; finally compute the total magnitude in all the other bands
; using the 4" and 6" aperture colors
    good = where((ndwfs.itot gt 0.0) and (ndwfs.bwcolor_4 gt -900.0) and (ndwfs.bwcolor_6 gt -900.0))
    ndwfs[good].bwtot_4     = ndwfs[good].itot + ndwfs[good].bwcolor_4
    ndwfs[good].bwtot_4_err = ndwfs[good].bwcolor_4_err
    ndwfs[good].bwtot_6     = ndwfs[good].itot + ndwfs[good].bwcolor_6
    ndwfs[good].bwtot_6_err = ndwfs[good].bwcolor_6_err

    good = where((ndwfs.itot gt 0.0) and (ndwfs.rcolor_4 gt -900.0) and (ndwfs.rcolor_6 gt -900.0))
    ndwfs[good].rtot_4     = ndwfs[good].itot + ndwfs[good].rcolor_4
    ndwfs[good].rtot_4_err = ndwfs[good].rcolor_4_err
    ndwfs[good].rtot_6     = ndwfs[good].itot + ndwfs[good].rcolor_6
    ndwfs[good].rtot_6_err = ndwfs[good].rcolor_6_err

    good = where((ndwfs.itot gt 0.0) and (ndwfs.kcolor_4 gt -900.0) and (ndwfs.kcolor_6 gt -900.0))
    ndwfs[good].ktot_4     = ndwfs[good].itot + ndwfs[good].kcolor_4
    ndwfs[good].ktot_4_err = ndwfs[good].kcolor_4_err
    ndwfs[good].ktot_6     = ndwfs[good].itot + ndwfs[good].kcolor_6
    ndwfs[good].ktot_6_err = ndwfs[good].kcolor_6_err

; for testing the K-band magnitudes
    if keyword_set(allndwfs) then begin
       good = where($
         (kband.k_mag_aper_03 gt 0.0) and (kband.k_mag_aper_03 lt 90.0) and $
         (kband.k_mag_aper_05 gt 0.0) and (kband.k_mag_aper_05 lt 90.0) and $
         (rband.r_mag_aper_03 gt 0.0) and (rband.r_mag_aper_03 lt 90.0) and $
         (rband.r_mag_aper_05 gt 0.0) and (rband.r_mag_aper_05 lt 90.0))
       ndwfs[good].kcolor_4_r = kband[good].k_mag_aper_03-rband[good].r_mag_aper_03 + (zpk-zpr)
       ndwfs[good].kcolor_6_r = kband[good].k_mag_aper_05-rband[good].r_mag_aper_05 + (zpk-zpr)
    endif else begin
       good = where($
         (kband.k_mag_aper3 gt 0.0) and (kband.k_mag_aper3 lt 90.0) and $
         (kband.k_mag_aper5 gt 0.0) and (kband.k_mag_aper5 lt 90.0) and $
         (rband.r_mag_aper3 gt 0.0) and (rband.r_mag_aper3 lt 90.0) and $
         (rband.r_mag_aper5 gt 0.0) and (rband.r_mag_aper5 lt 90.0))
       ndwfs[good].kcolor_4_r = kband[good].k_mag_aper3-rband[good].r_mag_aper3 + (zpk-zpr)
       ndwfs[good].kcolor_6_r = kband[good].k_mag_aper5-rband[good].r_mag_aper5 + (zpk-zpr)
    endelse

    good = where((ndwfs.rmag_auto gt 0.0) and (ndwfs.kcolor_4_r gt -900.0) and $
      (ndwfs.kcolor_6_r gt -900.0))
    ndwfs[good].ktot_r_4 = ndwfs[good].rmag_auto + ndwfs[good].kcolor_4_r
    ndwfs[good].ktot_r_6 = ndwfs[good].rmag_auto + ndwfs[good].kcolor_6_r
    
return, ndwfs
end

