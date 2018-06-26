;
; NAME:
;   GAMA_KCORRECT
;
; PURPOSE:
;   Compute K-corrections and rest-frame colors and absolute magnitudes for
;   GAMA/DR3 galaxies.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   clobber - overwrite an existing catalog
;
; OUTPUTS:
;   result - 
;
; COMMENTS:
;   No paths; filenames are hard-coded.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2018 Jun 25, Siena
;-

pro gama_kcorrect, kcorr, clobber=clobber

    ugriz_band_shift = 0.1 ; z=0.1
    vname = 'default.nolines' ; 'default'

    photfile = 'decals-dr5.0-GAMA-DR3-SpecObj-trim.fits'
    zcatfile = 'GAMA-DR3-SpecObj-trim.fits'

    phot = mrdfits(photfile,1)
    zcat = mrdfits(zcatfile,1)
    ngal = n_elements(zcat)

; only compute k-corrections for objects with good redshifts and photometry 
    good = where((zcat.nq ge 3) and (zcat.z ge 0.01) and (zcat.z le 1) and $
      (phot.flux_g gt 0) and (phot.flux_ivar_g gt 0) and $
      (phot.flux_r gt 0) and (phot.flux_ivar_r gt 0) and $
      (phot.flux_z gt 0) and (phot.flux_ivar_z gt 0),ngood)
;   good = good[0:100] & ngood = n_elements(good)
    
; convert the observed photometry to maggies
    splog, 'Converting observed photometry to maggies'
    obsfilters = legacysurvey_filterlist()
    nobsfilter = n_elements(obsfilters)
    
    obsmaggies = fltarr(nobsfilter,ngood)
    obsmaggies_ivar = fltarr(nobsfilter,ngood)

    obsmaggies[0,*] = 1e-9 * phot[good].flux_g / phot[good].mw_transmission_g
    obsmaggies[1,*] = 1e-9 * phot[good].flux_r / phot[good].mw_transmission_r
    obsmaggies[2,*] = 1e-9 * phot[good].flux_z / phot[good].mw_transmission_z
    
    obsmaggies_ivar[0,*] = phot[good].flux_ivar_g * phot[good].mw_transmission_g^2 * 1e18
    obsmaggies_ivar[1,*] = phot[good].flux_ivar_r * phot[good].mw_transmission_r^2 * 1e18
    obsmaggies_ivar[2,*] = phot[good].flux_ivar_z * phot[good].mw_transmission_z^2 * 1e18

; compute k-corrections    
    splog, 'Computing ugriz k-corrections'
    ugriz_kcorrect = im_kcorrect(zcat[good].z,obsmaggies,obsmaggies_ivar,$
      obsfilters,sdss_filterlist(),band_shift=ugriz_band_shift,chi2=chi2,$
      coeffs=coeffs,bestmaggies=bestmaggies,vname=vname,mass=mass,/silent,$
      absmag=ugriz_absmag,ivarabsmag=ugriz_absmag_ivar,clineflux=cflux)

    kcorr = {$
      cataid:                      -999L, $
      z:                          -999.0, $
      maggies:        fltarr(nobsfilter), $
      ivarmaggies:    fltarr(nobsfilter), $
      bestmaggies:    fltarr(nobsfilter), $
      mass:                       -999.0, $
;     intsfh:                     -999.0, $
      coeffs:                  fltarr(5), $
      chi2:                       -999.0, $
      cflux_4861:                 -999.0, $
      cflux_6563:                 -999.0, $
      ugriz_absmag_01:      fltarr(5)-999.0, $
      ugriz_absmag_01_ivar: fltarr(5)-999.0, $
      ugriz_kcorrect_01:    fltarr(5)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.cataid            = zcat.cataid
    kcorr.z                 = zcat.z
    kcorr[good].maggies     = obsmaggies
    kcorr[good].ivarmaggies = obsmaggies_ivar
    kcorr[good].bestmaggies = bestmaggies
    kcorr[good].mass        = alog10(mass)   ; Chabrier IMF
    kcorr[good].coeffs      = coeffs
    kcorr[good].chi2        = chi2
    kcorr[good].ugriz_absmag_01      = ugriz_absmag
    kcorr[good].ugriz_absmag_01_ivar = ugriz_absmag_ivar
    kcorr[good].ugriz_kcorrect_01    = ugriz_kcorrect

    kcorr[good].cflux_4861 = reform(cflux[1,*])
    kcorr[good].cflux_6563 = reform(cflux[4,*])
    
    outfile_kcorr = photfile.replace('.fits','-kcorr.fits')
    im_mwrfits, kcorr, outfile_kcorr, clobber=clobber
       
return
end
