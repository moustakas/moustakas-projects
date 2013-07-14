;
; NAME:
;       DEEP2_KCORRECT
;
; PURPOSE:
;       Compute K-corrections and rest-frame luminosities (and stellar
;       masses) for DEEP2 using K-correct.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       nokband - do not use the K-band photometry
;       write - write out the data structure
;
; OUTPUTS:
;       alldeep2phot - photometry data structure
;       allresult   - stellar masses and K-corrections
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2008 Mar 25, NYU - stolen from AGES_KCORRECT_ALL
;       jm08sep04nyu - updated to use IM_KCORRECT()
;       jm13jul13siena - updated to most recent version of
;         IM_KCORRECT() 
;-

pro deep2_kcorrect, zcat, result, clobber=clobber

; define the cosmology and some constants

    version = deep2_version(/kcorr)
    analysis_path = deep2_path(/analysis)

    omega0 = 0.3 & omegal0 = 0.7 & h100 = 0.7
    
    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0

;   vname = 'default'
    vname = 'default.nolines' ; 'default'

; read the merged photometric catalog; only compute k-corrections for
; objects with well-measured redshifts
    if (n_elements(zcat) eq 0L) then zcat = read_deep2_zcat()
    ngalaxy = n_elements(zcat)

; convert the observed photometry to maggies
    splog, 'Converting observed photometry to maggies'
    obsfilters = deep2_filterlist()
    nobsfilter = n_elements(obsfilters)
    deep2_to_maggies, zcat, obsmaggies, obsmaggies_ivar

; compute k-corrections    
    
    splog, 'Computing ugriz k-corrections'
    ugriz_kcorrect = im_kcorrect(zcat.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,sdss_filterlist(),band_shift=ugriz_band_shift,chi2=chi2,$
      coeffs=coeffs,rmaggies=rmaggies,vname=vname,mass=mass,$
      absmag=ugriz_absmag,ivarabsmag=ugriz_absmag_ivar,$
      clineflux=cflux,omega0=omega0,omegal0=omegal0,/silent)

    splog, 'Computing UBVRI k-corrections'
    ubvri_kcorrect = im_kcorrect(zcat.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,band_shift=ubvri_band_shift,/vega,$ ; note /VEGA
      vname=vname,absmag=ubvri_absmag,ivarabsmag=ubvri_absmag_ivar,$
      omega0=omega0,omegal0=omegal0,/silent)

    result_template = {$
;     z:                          -999.0, $
      abmaggies:      fltarr(nobsfilter), $
      abmaggies_ivar: fltarr(nobsfilter), $
      mass:                       -999.0, $
      intsfh:                     -999.0, $
      coeffs:                  fltarr(5), $
      chi2:                       -999.0, $
      cflux_3727:                 -999.0, $
      cflux_4861:                 -999.0, $
      cflux_4959:                 -999.0, $
      cflux_5007:                 -999.0, $
      cflux_6563:                 -999.0, $
      ugriz_absmag:      fltarr(5)-999.0, $
      ugriz_absmag_ivar: fltarr(5)-999.0, $
      ugriz_kcorrect:    fltarr(5)-999.0, $
      ubvri_absmag:      fltarr(5)-999.0, $
      ubvri_absmag_ivar: fltarr(5)-999.0, $
      ubvri_kcorrect:    fltarr(5)-999.0}
    result = replicate(result_template,ngalaxy)
    result = struct_addtags(zcat,temporary(result))

;   result.z                 = zcat.z
    result.abmaggies         = obsmaggies
    result.abmaggies_ivar    = obsmaggies_ivar
    result.mass              = alog10(mass)   ; Chabrier IMF
;   result.intsfh            = alog10(intsfh) ; Chabrier IMF
    result.coeffs            = coeffs
    result.chi2              = chi2
    result.ugriz_absmag      = ugriz_absmag
    result.ugriz_absmag_ivar = ugriz_absmag_ivar
    result.ugriz_kcorrect    = ugriz_kcorrect
    result.ubvri_absmag      = ubvri_absmag
    result.ubvri_absmag_ivar = ubvri_absmag_ivar
    result.ubvri_kcorrect    = ubvri_kcorrect

    result.cflux_3727 = reform(cflux[0,*])
    result.cflux_4861 = reform(cflux[1,*])
    result.cflux_4959 = reform(cflux[2,*])
    result.cflux_5007 = reform(cflux[3,*])
    result.cflux_6563 = reform(cflux[4,*])
    
    outfile_kcorr = analysis_path+'deep2_kcorr_'+version+'.fits'
    splog, 'Writing '+outfile_kcorr
    im_mwrfits, result, outfile_kcorr, clobber=clobber
       
return
end

;;; compute stellar masses according to Bell & de Jong; need to adopt
;;; the solar magnitudes in the Bell & de Jong paper!
;;
;;    bell03_johnson = read_03bell(/johnson)
;;
;;;   absmsun_r = 4.46            
;;    absmsun_r = k_solar_magnitudes(filterlist='bessell_R.par')
;;
;;    bellcoeff = [bell03_johnson[1].ar,bell03_johnson[1].br] ; [bell01[1].ar,bell01[1].br]
;;    color = (result.ubvri_absmag[1]-result.ubvri_absmag[3])
;;    result.ml_br_r = (bellcoeff[0] + bellcoeff[1]*color)*(h100/h100_bell03) ; h=1-->h=0.7
;;    result.mass_br_r = result.ml_br_r + (-0.4*(result.ubvri_absmag[3]-absmsun_r)) ; h=0.7

