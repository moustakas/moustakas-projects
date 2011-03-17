pro init_sfrf_results, nzbins, sfrf_fit=sfrf_fit, sfrf_data=sfrf_data, $
  double_pl=double_pl
; jm10mar19ucsd - initialize the SFRF results (the SFRFs themselves
;   and the best-fitting parameters structure)

; best-fit parameters    
    sfrf_fit = {$
      pl:         keyword_set(double_pl), $
      minsfr:      0.0,$
      phistar:     0.0,$
      sfrstar:     0.0,$
      alpha:       0.0,$
      phiplus:     0.0,$
      alphaplus:   0.0,$

      phistar_err:   0.0,$
      sfrstar_err:   0.0,$
      alpha_err:     0.0,$
      phiplus_err:   0.0,$
      alphaplus_err: 0.0,$

      phistar_cv_err:   0.0,$ ; cosmic variance
      sfrstar_cv_err:     0.0,$ 
      alpha_cv_err:     0.0,$
      phiplus_cv_err:   0.0,$
      alphaplus_cv_err: 0.0,$

      phistar_model_err:   0.0,$ ; various stellar sfres
      sfrstar_model_err:     0.0,$ 
      alpha_model_err:     0.0,$
      phiplus_model_err:   0.0,$
      alphaplus_model_err: 0.0}
    sfrf_fit = replicate(sfrf_fit,nzbins)

; MF data
    null = fltarr(50)-999.0
    sfrf_data = {$
      ngal:           0L,$
      nbins:           0,$
      fullbin: fix(null),$
      number:  fix(null),$
      sfr:         null,$
      phi:          null,$
      phierr:       null,$
      phierr_cv:    null,$ ; cosmic variance
      phierr_model: null}  ; assuming various stellar sfres
    sfrf_data = replicate(sfrf_data,nzbins)

return
end
