pro redmapper_kcorrect, clobber=clobber
; jm13mar28siena - compute K-corrections for the REDMAPPER project

;   kcorr = mrdfits(redmapperpath+'redmapper_kcorrect.fits.gz',1)
;   kcorrect_qaplot, kcorr, psfile=redmapperpath+'qaplot_redmapper_kcorrect.ps', $
;     in_filterlist=redmapper_filterlist(), /clobber, vname='default.nolines'

    redmapperpath = redmapper_path(/catalogs,version=ver)

    h100 = 0.7
    vname = 'default'
;   vname = 'default.nolines'

    kcorrfile = redmapperpath+'redmapper_'+ver+'_kcorrect.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    filters = redmapper_filterlist()
    nfilt = n_elements(filters)

    cat = mrdfits(redmapperpath+'redmapper_'+ver+'_photometry.fits.gz',1)
    ngal = n_elements(cat)

    kcorr = {$
      zobj:                          -999.0,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      mass:                          -999.0,$
      coeffs:                     fltarr(5),$
      chi2:                          -999.0,$
      uvflux:                     fltarr(2),$

      mobs_fnuvugrizjhk:              fltarr(10)-999.0,$
      fnuvugrizjhk_absmag_00:         fltarr(10)-999.0,$
      fnuvugrizjhk_absmag_ivar_00:    fltarr(10)-999.0,$
      fnuvugrizjhk_kcorrect_00:       fltarr(10)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.zobj = cat.z
    kcorr.maggies = cat.maggies
    kcorr.ivarmaggies = cat.ivarmaggies
    
; compute k-corrections
    out_filters = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]
    
    splog, 'Computing [FN]UV-ugriz-JHK K-corrections'
    kcorrect_00 = im_kcorrect(cat.z,cat.maggies,cat.ivarmaggies,$
      filters,out_filters,band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag_00,$
      ivarabsmag=absmag_ivar_00,uvflux=uvflux,$;clineflux=cflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.0

; get the apparent and absolute magnitudes in ugrizJHKs
    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_reconstruct_maggies, coeffs, cat.z, appmaggies, vmatrix=vmatrix, $
      lambda=lambda, filterlist=out_filters
    kcorr.mobs_fnuvugrizjhk = reform(-2.5*alog10(appmaggies))
    
    kcorr.bestmaggies = bestmaggies
    kcorr.mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.coeffs = coeffs
    kcorr.chi2 = chi2
    kcorr.uvflux = uvflux
    kcorr.fnuvugrizjhk_absmag_00      = absmag_00
    kcorr.fnuvugrizjhk_absmag_ivar_00 = absmag_ivar_00
    kcorr.fnuvugrizjhk_kcorrect_00    = kcorrect_00

    im_mwrfits, kcorr, kcorrfile, /clobber

return
end
