pro init_mf_results, nzbins, mf_fit=mf_fit, mf_data=mf_data, $
  double_schechter=double_schechter
; jm10feb15ucsd - initialize the MF results (the MFs themselves and
;   the best-fitting parameters structure)

; best-fit parameters    
    mf_fit = {$
      double:     keyword_set(double_schechter), $
      minmass:     0.0,$
      phistar:     0.0,$
      mstar:       0.0,$
      alpha:       0.0,$
      phiplus:     0.0,$
      alphaplus:   0.0,$

      phistar_err:   0.0,$
      mstar_err:     0.0,$
      alpha_err:     0.0,$
      phiplus_err:   0.0,$
      alphaplus_err: 0.0,$

      phistar_cv_err:   0.0,$ ; cosmic variance
      mstar_cv_err:     0.0,$ 
      alpha_cv_err:     0.0,$
      phiplus_cv_err:   0.0,$
      alphaplus_cv_err: 0.0,$

      phistar_model_err:   0.0,$ ; various stellar masses
      mstar_model_err:     0.0,$ 
      alpha_model_err:     0.0,$
      phiplus_model_err:   0.0,$
      alphaplus_model_err: 0.0}
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
      phierr_cv:    null,$ ; cosmic variance
      phierr_model: null}  ; assuming various stellar masses
    mf_data = replicate(mf_data,nzbins)

return
end
