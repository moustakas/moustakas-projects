; ###########################################################################
; OBSOLETE! - SEE deep2_kcorrect
; ###########################################################################
stop
pro deep2_kcorrect_all, zcat, result, write=write
; jm08mar25nyu - stolen from AGES_KCORRECT_ALL
    
; define the cosmology and some constants

    version = deep2_version(/kcorr)
    analysis_path = deep2_path(/analysis)

    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()
    
    h100_kcorrect = 1.0
    h100_bell01 = 0.71
    h100_bell03 = 1.0

    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0

    vname = 'default.nolines' ; 'default'

; read the DEEP2 redshift catalog; THIS SHOULD MATCH DEEP2_SPECFIT!! 

    if (n_elements(zcat) eq 0L) then zcat = read_deep2_zcat(/good,dr=dr) 
    ngalaxy = n_elements(zcat)

; convert the observed photometry to maggies and store

    splog, 'Converting observed photometry to maggies.'
    obsfilters = ['deep_'+['B','R','I']]+'.par'
    nobsfilter = n_elements(obsfilters)
    deep_to_maggies, zcat, obsmaggies, obsmaggies_ivar

; compute k-corrections    

    splog, 'Computing ugriz k-corrections.'
    ugriz_kcorrect = deep2_kcorrect(zcat.z,nmgy=obsmaggies*1E9,ivar=obsmaggies_ivar*1E-18,$
      /closest,band_shift=ugriz_band_shift,/sdss,$ ; note /SDSS
      chi2=ugriz_chi2,coeffs=ugriz_coeffs,rmaggies=ugriz_rmaggies,omaggies=ugriz_omaggies,$
      oivar=ugriz_oivar,vname=vname,mass=ugriz_mass,absmag=ugriz_absmag,$
      amivar=ugriz_amivar,omega0=omega0,omegal0=omegal,intsfh=ugriz_intsfh,/silent)

    splog, 'Computing UBVRI k-corrections.'
    ubvri_kcorrect = deep2_kcorrect(zcat.z,nmgy=obsmaggies*1E9,ivar=obsmaggies_ivar*1E-18,$
      /closest,band_shift=ubvri_band_shift,/vega,$ ; note /VEGA
      chi2=ubvri_chi2,coeffs=ubvri_coeffs,rmaggies=ubvri_rmaggies,omaggies=ubvri_omaggies,$
      oivar=ubvri_oivar,vname=vname,mass=ubvri_mass,absmag=ubvri_absmag,$
      amivar=ubvri_amivar,omega0=omega0,omegal0=omegal,intsfh=ubvri_intsfh,/silent)

    result_template = {$
;     z:                             -999.0, $
      abmaggies:         fltarr(nobsfilter), $
      abmaggies_ivar:    fltarr(nobsfilter), $
      mass:                          -999.0, $
      intsfh:                        -999.0, $
      coeffs:                     fltarr(5), $
      cflux_3727:                    -999.0, $
      cflux_4861:                    -999.0, $
      cflux_4959:                    -999.0, $
      cflux_5007:                    -999.0, $
      cflux_6563:                    -999.0, $
      ugriz_absmag:         fltarr(5)-999.0, $
      ugriz_absmag_ivar:    fltarr(5)-999.0, $
      ugriz_kcorrect:       fltarr(5)-999.0, $
      ubvri_absmag:         fltarr(5)-999.0, $
      ubvri_absmag_ivar:    fltarr(5)-999.0, $
      ubvri_kcorrect:       fltarr(5)-999.0}
;     ml_br_r:                       -999.0, $
;     mass_br_r:                     -999.0}
    result = replicate(result_template,ngalaxy)
    result = struct_addtags(zcat,temporary(result))

;   result.z                 = zcat.z
    result.abmaggies         = ugriz_omaggies
    result.abmaggies_ivar    = ugriz_oivar
    result.mass              = alog10(ugriz_mass) - 2.0*alog10(h100/h100_kcorrect)   ; h=1-->h=0.7 (Chabrier IMF)
    result.intsfh            = alog10(ugriz_intsfh) - 2.0*alog10(h100/h100_kcorrect) ; h=1-->h=0.7 (Chabrier IMF)
    result.coeffs            = ugriz_coeffs
    result.ugriz_absmag      = ugriz_absmag + 5.0*alog10(h100/h100_kcorrect) ; h=1-->h=0.7
    result.ugriz_absmag_ivar = ugriz_amivar
    result.ugriz_kcorrect    = ugriz_kcorrect
    result.ubvri_absmag      = ubvri_absmag + 5.0*alog10(h100/h100_kcorrect) ; h=1-->h=0.7
    result.ubvri_absmag_ivar = ubvri_amivar
    result.ubvri_kcorrect    = ubvri_kcorrect

    splog, 'Computing emission-line continuum fluxes.'
    k_load_vmatrix, vmatrix, lambda, vname=vname
    restwave = k_lambda_to_centers(lambda)

    linewave = [3727.420,4861.325,4958.911,5006.843,6562.800]
    lineindx = findex(restwave,linewave)
    
    cflux = fltarr(5,ngalaxy)
    for ii = 0L, ngalaxy-1L do cflux[*,ii] = $
      interpolate(reform(vmatrix#ugriz_coeffs[*,ii]),lineindx) ; [erg/s/cm2/A]
    
    result.cflux_3727 = reform(cflux[0,*])
    result.cflux_4861 = reform(cflux[1,*])
    result.cflux_4959 = reform(cflux[2,*])
    result.cflux_5007 = reform(cflux[3,*])
    result.cflux_6563 = reform(cflux[4,*])
    
; compute stellar masses according to Bell & de Jong; need to adopt
; the solar magnitudes in the Bell & de Jong paper!
;
;   bell03_johnson = read_03bell(/johnson)
;
;   absmsun_r = 4.46            
;   absmsun_r = k_solar_magnitudes(filterlist='bessell_R.par')
;
;   bellcoeff = [bell03_johnson[1].ar,bell03_johnson[1].br] ; [bell01[1].ar,bell01[1].br]
;   color = (result.ubvri_absmag[1]-result.ubvri_absmag[3])
;   result.ml_br_r = bellcoeff[0] + bellcoeff[1]*color
;   result.mass_br_r = result.ml_br_r + $
;     (-0.4*(result.ubvri_absmag[3]-absmsun_r)) - 2.0*alog10(h100/h100_bell03) ; h=1-->h=0.7

    if keyword_set(write) then begin
       outfile_kcorr = analysis_path+'deep2_kcorr_'+version+'.fits'
       splog, 'Writing '+outfile_kcorr
       mwrfits, result, outfile_kcorr, /create
       spawn, 'gzip -f '+outfile_kcorr, /sh
    endif
       
return
end
