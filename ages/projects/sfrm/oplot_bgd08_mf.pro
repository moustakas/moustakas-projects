pro oplot_bgd08_mf, maxis, params=params, mf=mf, log=log, $
  salpeter=salpeter, _extra=extra
; jm10feb15ucsd - overplot the Baldry, Glazebrook, & Driver (2008)
; stellar mass function

    if (n_elements(maxis) eq 0) then maxis = im_array(7.5,12.5,0.01)

    params = {$
      phistar:   4.26E-3,$
      mstar:      10.648,$
      alpha:       -0.46,$
      phiplus:   0.58E-3,$
      alphaplus:   -1.58}

    if keyword_set(salpeter) then offset = +0.14 else offset = 0.0
;   if keyword_set(salpeter) then offset = -alog10(0.7) else offset = 0.0

;   params = [4.26E-3,10.0^10.648,-0.46,0.58E-3,-1.58]
    mf = mf_schechter_plus(maxis,params)
    if keyword_set(log) then mf = alog10(mf)
    djs_oplot, maxis+offset, mf, _extra=extra

return
end
    
