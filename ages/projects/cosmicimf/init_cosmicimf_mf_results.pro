pro init_cosmicimf_mf_results, nzbins, mf_fit=mf_fit, mf_data=mf_data
; jm10mar23ucsd - initialize the MF results (the MFs themselves and
;   the best-fitting parameters structure)

; best-fit parameters    
    mf_fit = {$
      chi2_dof:    0.0,$ ; chi^2_nu
      
      phistar:     0.0,$
      mstar:       0.0,$
      alpha:       0.0,$

      phistar_err:   0.0,$
      mstar_err:     0.0,$
      alpha_err:     0.0,$

      phistar_cv_err:   0.0,$ ; cosmic variance
      mstar_cv_err:     0.0,$ 
      alpha_cv_err:     0.0,$

      phistar_lss_err:   0.0,$ ; LSS on even larger scales

      rho:           0.0,$  ; stellar mass density [M_sun Mpc^3]
      rho_err:       0.0,$
      rho_cv_err:    0.0,$
      rho_lss_err:   0.0}
    mf_fit = replicate(mf_fit,nzbins)

; MF data
    null = fltarr(50)-999.0
    mf_data = {$
      ngal:           0L,$
      nbins:           0,$
      fullbin: fix(null),$
      number:  fix(null),$
      mass:         null,$
      phi:          null,$
      phierr:       null,$
      phierr_cv:    null}  ; cosmic variance
;     phierr_model: null}  ; assuming various stellar masses
    mf_data = replicate(mf_data,nzbins)

return
end
