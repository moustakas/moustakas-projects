pro init_lf24_results, nzbins, nimf, lf24_fit=lf24_fit, lf24_data=lf24_data
; jm10mar25ucsd - initialize the L(24) LF results (the LFs themselves
;   and the best-fitting parameters structure)

; best-fit parameters    
    lf24_fit = {$
      chi2_dof:    0.0,$ ; chi^2_nu

      phistar:     0.0,$
      lstar:       0.0,$
      alpha:       0.0,$
      beta:        0.0,$

      phistar_err:   0.0,$
      lstar_err:     0.0,$
      alpha_err:     0.0,$
      beta_err:      0.0,$

      lstar_cv_err:  0.0,$

      rhosfr:           fltarr(nimf),$ ; SFR density for each IMF [M_sun/yr/Mpc^3]
      rhosfr_err:       fltarr(nimf),$
      rhosfr_lss_err:   fltarr(nimf),$
      rhosfr_cv_err:    fltarr(nimf),$
      
      rhosfr_local:     fltarr(nimf),$ ; local SFR density for each IMF [M_sun/yr/Mpc^3]
      rhosfr_local_err: fltarr(nimf)}
    lf24_fit = replicate(lf24_fit,nzbins)

; LF data
    null = fltarr(50)-999.0
    lf24_data = {$
      ngal:           0L,$
      nbins:           0,$
      fullbin: fix(null),$
      number:  fix(null),$
      l24:          null,$
      phi:          null,$
      phierr:       null,$
      phierr_cv:    null} ; cosmic variance
    lf24_data = replicate(lf24_data,nzbins)

return
end
