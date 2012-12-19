function test_isedfit_kcorrect, in_redshift, in_maggies, in_ivarmaggies, $
  filterlist=filterlist
; jm10aug30ucsd - simple wrapper script to compute a general set of
; K-corrections for the PRIMUS/MF project (called by MF_ISEDFIT)

    h100 = 0.7
    vname = 'default'
;   vname = 'default.nolines'
    
    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    kcorr = {$
      k_zobj:                          -999.0, $
      k_maggies:                fltarr(nfilt), $
      k_ivarmaggies:            fltarr(nfilt), $
      k_bestmaggies:            fltarr(nfilt), $
      k_mass:                          -999.0, $
      k_coeffs:                     fltarr(5), $
      k_chi2:                          -999.0, $
      k_uvflux:                      fltarr(2),$

      k_galex_absmag_01:         fltarr(2)-999.0,$
      k_galex_absmag_ivar_01:    fltarr(2)-999.0,$
      k_galex_kcorrect_01:       fltarr(2)-999.0,$

      k_ugriz_absmag_01:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_01:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_01:       fltarr(5)-999.0,$

      k_ubvri_absmag_00:         fltarr(5)-999.0,$
      k_ubvri_absmag_ivar_00:    fltarr(5)-999.0,$
      k_ubvri_kcorrect_00:       fltarr(5)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = in_redshift
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
    
; compute k-corrections    
    splog, 'Computing ugriz K-corrections'
    ugriz_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[galex_filterlist(),sdss_filterlist()],band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag_01,$
      ivarabsmag=ugriz_absmag_ivar_01,/silent,vname=vname,h100=h100,$
      uvflux=uvflux)                         ; AB, band_shift=0.1

    splog, 'Computing ubvri K-corrections'
    ubvri_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,bessell_filterlist(),band_shift=0.0,absmag=ubvri_absmag_00,$
      ivarabsmag=ubvri_absmag_ivar_00,/silent,vname=vname,h100=h100) ; AB, band_shift=0.0

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_uvflux = uvflux

    kcorr.k_galex_absmag_01      = ugriz_absmag_01[0:1,*]
    kcorr.k_galex_absmag_ivar_01 = ugriz_absmag_ivar_01[0:1,*]
    kcorr.k_galex_kcorrect_01    = ugriz_kcorrect_01[0:1,*]

    kcorr.k_ugriz_absmag_01      = ugriz_absmag_01[2:6,*]
    kcorr.k_ugriz_absmag_ivar_01 = ugriz_absmag_ivar_01[2:6,*]
    kcorr.k_ugriz_kcorrect_01    = ugriz_kcorrect_01[2:6,*]

    kcorr.k_ubvri_absmag_00      = ubvri_absmag_00
    kcorr.k_ubvri_absmag_ivar_00 = ubvri_absmag_ivar_00
    kcorr.k_ubvri_kcorrect_00    = ubvri_kcorrect_00

return, kcorr
end

