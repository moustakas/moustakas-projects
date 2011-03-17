function fit_irx_a1500, irx, a1500, model_irx=model_irx, $
  model_a1500=model_a1500, debug=debug
; fit the A(1500) vs IRX relation using eq (2) in Gordon+00

    xx = irx
    xx1 = interpol(xx,a1500,1.0)
    xx2 = interpol(xx,a1500,1.75)
    
    lo = where(xx lt xx1)
    coeff_lo = poly_fit(xx[lo],a1500[lo],3,yfit=yfit)
    hi = where((xx gt xx2) and (xx lt 1E3))
    coeff_hi = poly_fit(alog(xx[hi]),a1500[hi],2,yfit=yfit)
;   djs_plot, xx[hi], a1500[hi], ps=6

    coeff = [xx1,xx2,reform(coeff_lo),reform(coeff_hi)] ; final coefficients [9]
    
; interpolate onto a desired output grid    
    if (n_elements(model_irx) ne 0) then begin
       model_a1500_lo = poly(model_irx,coeff_lo)
       model_a1500_hi = poly(alog(model_irx),coeff_hi)
       model_ww = (xx2-model_irx)/(xx2-xx1)
       model_a1500_mid = model_ww*model_a1500_lo+(1-model_ww)*model_a1500_hi

       mlo = where(model_irx lt xx1)
       mmid = where((model_irx gt xx1) and (model_irx lt xx2))
       mhi = where(model_irx gt xx2)
       model_a1500 = model_irx*0.0
       model_a1500[mlo] = poly(model_irx[mlo],coeff_lo)
       model_a1500[mhi] = poly(alog(model_irx[mhi]),coeff_hi)

       model_ww_mid = (xx2-model_irx[mmid])/(xx2-xx1)
       model_ax_mid = poly(model_irx[mmid],coeff_lo)
       model_bx_mid = poly(alog(model_irx[mmid]),coeff_hi)
       
       model_a1500[mmid] = model_ww_mid*model_ax_mid+(1-model_ww_mid)*model_bx_mid

       if keyword_set(debug) then begin
          djs_plot, irx, a1500, psym=6, xrange=[0.1,1E3], yr=[0,8], /xlog
          djs_oplot, model_irx, model_a1500_lo, color='yellow'
          djs_oplot, model_irx, model_a1500_mid, color='green'
          djs_oplot, model_irx, model_a1500_hi, color='red'
          djs_oplot, model_irx, model_a1500, color='cyan', thick=3
       endif
    endif
    
return, coeff
end

