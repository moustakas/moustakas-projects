pro redmapper_kcorrect, clobber=clobber
; jm13mar28siena - compute K-corrections for the REDMAPPER project

;   kcorr = mrdfits(redmapperpath+'redmapper_kcorrect.fits.gz',1)
;   kcorrect_qaplot, kcorr, psfile=redmapperpath+'qaplot_redmapper_kcorrect.ps', $
;     in_filterlist=redmapper_filterlist(), /clobber, vname='default.nolines'

    catpath = redmapper_path(version=ver,/catalogs)

    h100 = 0.7
    vname = 'default'
;   vname = 'default.nolines'

    kcorrfile = catpath+'redmapper_'+ver+'_kcorr.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    filters = redmapper_filterlist()
    nfilt = n_elements(filters)

    cat = mrdfits(catpath+'redmapper_'+ver+'_phot.fits.gz',1);,range=[0,10])
    ngal = n_elements(cat)

    kcorr = {$
      z:                             -999.0,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      kcorr_mstar:                   -999.0,$
      kcorr_coeffs:               fltarr(5),$
      kcorr_chi2:                    -999.0,$
      kcorr_uvflux:               fltarr(2),$

;     mobs_fnuvugrizjhk:           fltarr(10)-999.0,$
      fnuv_absmag_00:               fltarr(2)-999.0,$
      fnuv_absmag_ivar_00:          fltarr(2)-999.0,$
      ugriz_absmag_00:              fltarr(5)-999.0,$
      ugriz_absmag_ivar_00:         fltarr(5)-999.0,$
      ubvri_absmag_00:              fltarr(5)-999.0,$
      ubvri_absmag_ivar_00:         fltarr(5)-999.0,$
      jhk_absmag_00:                fltarr(3)-999.0,$
      jhk_absmag_ivar_00:           fltarr(3)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.z = cat.z
    kcorr.maggies = cat.maggies
    kcorr.ivarmaggies = cat.ivarmaggies
    
; compute k-corrections
    out_filters = [galex_filterlist(),sdss_filterlist(),twomass_filterlist(),bessell_filterlist()]
    
    splog, 'Computing [FN]UV-ugriz-JHK K-corrections'
    kcorrect_00 = im_kcorrect(cat.z,cat.maggies,cat.ivarmaggies,$
      filters,out_filters,band_shift=0.0,chi2=chi2,mass=mstar,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag_00,$
      ivarabsmag=absmag_ivar_00,uvflux=uvflux,$;clineflux=cflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.0

;; get the apparent and absolute magnitudes in ugrizJHKs
;    k_load_vmatrix, vmatrix, lambda, vname=vname
;    k_reconstruct_maggies, coeffs, cat.z, appmaggies, vmatrix=vmatrix, $
;      lambda=lambda, filterlist=out_filters
;    kcorr.mobs_fnuvugrizjhk = reform(-2.5*alog10(appmaggies))
    
    kcorr.bestmaggies = bestmaggies
    kcorr.kcorr_mstar = alog10(mstar) ; h=0.7, Chabrier
    kcorr.kcorr_coeffs = coeffs
    kcorr.kcorr_chi2 = chi2
    kcorr.kcorr_uvflux = uvflux

    kcorr.fnuv_absmag_00       = reform(absmag_00[0:1,*])
    kcorr.fnuv_absmag_ivar_00  = reform(absmag_ivar_00[0:1,*])
    kcorr.ugriz_absmag_00      = reform(absmag_00[2:6,*])
    kcorr.ugriz_absmag_ivar_00 = reform(absmag_ivar_00[2:6,*])
    kcorr.ubvri_absmag_00      = reform(absmag_00[7:11,*])
    kcorr.ubvri_absmag_ivar_00 = reform(absmag_ivar_00[7:11,*])
    kcorr.jhk_absmag_00        = reform(absmag_00[12:14,*])
    kcorr.jhk_absmag_ivar_00   = reform(absmag_ivar_00[12:14,*])

    im_mwrfits, kcorr, kcorrfile, /clobber

stop    
    
return
end
