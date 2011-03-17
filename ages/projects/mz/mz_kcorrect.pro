function mz_kcorrect, in_redshift, in_maggies, in_ivarmaggies, $
  filterlist=filterlist, just_ugriz=just_ugriz
; simple wrapper script to compute a general set of K-corrections for
; the AGES/MZ project

    vname = mz_vname()
    h100 = mz_h100()
    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    
    kcorr = {$
      k_z:                             -999.0, $
      k_maggies:                fltarr(nfilt), $
      k_ivarmaggies:            fltarr(nfilt), $
      k_bestmaggies:            fltarr(nfilt), $
      k_mass:                          -999.0, $
      k_coeffs:                     fltarr(5), $
      k_chi2:                          -999.0, $
      k_cflux:                fltarr(5)-999.0, $
      k_uvflux:               fltarr(2)-999.0, $

      k_galex_absmag_00:         fltarr(2)-999.0,$
      k_galex_absmag_ivar_00:    fltarr(2)-999.0,$
      k_galex_kcorrect_00:       fltarr(2)-999.0,$
;     k_galex_absmag_01:         fltarr(2)-999.0,$
;     k_galex_absmag_ivar_01:    fltarr(2)-999.0,$
;     k_galex_kcorrect_01:       fltarr(2)-999.0,$

;     k_ugriz_absmag_00:         fltarr(5)-999.0,$
;     k_ugriz_absmag_ivar_00:    fltarr(5)-999.0,$
;     k_ugriz_kcorrect_00:       fltarr(5)-999.0,$
      k_ugriz_absmag_01:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_01:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_01:       fltarr(5)-999.0,$

      k_ubvrijhk_absmag_00:      fltarr(8)-999.0,$
      k_ubvrijhk_absmag_ivar_00: fltarr(8)-999.0,$
      k_ubvrijhk_kcorrect_00:    fltarr(8)-999.0}
;     k_ubvrijhk_absmag_01:      fltarr(8)-999.0,$
;     k_ubvrijhk_absmag_ivar_01: fltarr(8)-999.0,$
;     k_ubvrijhk_kcorrect_01:    fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_z = in_redshift
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
    
; compute k-corrections    
    splog, 'Computing ugriz K-corrections'
    ugriz_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag_01,$
      ivarabsmag=ugriz_absmag_ivar_01,clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.1

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_cflux = cflux
    kcorr.k_uvflux = uvflux
    kcorr.k_ugriz_absmag_01         = ugriz_absmag_01
    kcorr.k_ugriz_absmag_ivar_01    = ugriz_absmag_ivar_01
    kcorr.k_ugriz_kcorrect_01       = ugriz_kcorrect_01

    if (keyword_set(just_ugriz) eq 0) then begin
       splog, 'Computing NUV,FUV K-corrections'
       galex_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
         filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_00,$
         ivarabsmag=galex_absmag_ivar_00,/silent,vname=vname,h100=h100) ; AB, band_shift=0.0

       splog, 'Computing UBVRIJHKs K-corrections'
       kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
         filterlist,[bessell_filterlist(),twomass_filterlist()],band_shift=0.0,$
         absmag=absmag_00,ivarabsmag=absmag_ivar_00,/silent,vname=vname,h100=h100) ; AB, band_shift=0.0
       
       kcorr.k_galex_absmag_00         = galex_absmag_00
       kcorr.k_galex_absmag_ivar_00    = galex_absmag_ivar_00
       kcorr.k_galex_kcorrect_00       = galex_kcorrect_00

       kcorr.k_ubvrijhk_absmag_00      = absmag_00
       kcorr.k_ubvrijhk_absmag_ivar_00 = absmag_ivar_00
       kcorr.k_ubvrijhk_kcorrect_00    = kcorrect_00
    endif

return, kcorr
end

