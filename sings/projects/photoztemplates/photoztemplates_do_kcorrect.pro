function photoztemplates_do_kcorrect, in_redshift, in_maggies, in_ivarmaggies, $
  filterlist=filterlist, silent=silent
; jm11apr27ucsd - compute K-corrections

    h100 = 0.7
;   vname = 'default'
    vname = 'default.nolines'

    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)

    kcorr = {$
      k_zobj:                          -999.0,$
      k_maggies:                fltarr(nfilt),$
      k_ivarmaggies:            fltarr(nfilt),$
      k_bestmaggies:            fltarr(nfilt),$
      k_mass:                          -999.0,$
      k_coeffs:                     fltarr(5),$
      k_chi2:                          -999.0,$

      k_uvflux:               fltarr(2)-999.0, $

      k_galex_absmag:         fltarr(2)-999.0,$
      k_galex_absmag_ivar:    fltarr(2)-999.0,$
      k_galex_kcorrect:       fltarr(2)-999.0,$

      k_ugrizjhk_absmag:         fltarr(8)-999.0,$
      k_ugrizjhk_absmag_ivar:    fltarr(8)-999.0,$
      k_ugrizjhk_kcorrect:       fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = in_redshift
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
    
; compute k-corrections
    out_filterlist = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]
    
    if (keyword_set(silent) eq 0) then splog, 'Computing K-corrections'
    kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,out_filterlist,band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag,$
      ivarabsmag=absmag_ivar,uvflux=uvflux,/silent,vname=vname,h100=h100) ; AB, band_shift=0.0

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_uvflux = uvflux

    kcorr.k_galex_absmag         = absmag[0:1,*]
    kcorr.k_galex_absmag_ivar    = absmag_ivar[0:1,*]
    kcorr.k_galex_kcorrect       = kcorrect[0:1,*]

    kcorr.k_ugrizjhk_absmag         = absmag[2:9,*]
    kcorr.k_ugrizjhk_absmag_ivar    = absmag_ivar[2:9,*]
    kcorr.k_ugrizjhk_kcorrect       = kcorrect[2:9,*]

return, kcorr
end

