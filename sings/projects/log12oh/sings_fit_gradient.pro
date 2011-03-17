function sings_fit_gradient, x, xerr, y, yerr1, nmonte=nmonte, $
  xchar=xchar, xcent=xcent, noslopeconstrain=noslopeconstrain, $
  weightedfit=weightedfit, yminerr=yminerr, ohmodel=yymodel
; jm06mar29uofa - excised from SINGS_LOG12OH
; jm08feb11nyu - use iterative outlier rejection
; jm10jul06ucsd - minor tweaks
    
    npts = n_elements(x)
    if (npts eq 0L) then begin
       return, -1L
    endif
    
    if (n_elements(nmonte) eq 0) then nmonte = 100
    if (n_elements(xchar) eq 0) then xchar = 0.4
    if (n_elements(xcent) eq 0) then xcent = 0.0
    if (n_elements(yminerr) eq 0) then yminerr = 0.0

    neg = where((xerr lt 0.0) or (yerr1 lt 0.0),nneg)
    if (nneg ne 0L) then message, 'Negative errors!'

    yerr = sqrt(yerr1^2 + yminerr^2)
    rr25axis = findgen(501)/100.0

    res = {$
      slope_flag:      0,$ ; formal positive slope
      int:     [0.0,0.0],$
      slope:   [0.0,0.0],$
      ychar:   [0.0,0.0],$
      ycent:   [0.0,0.0],$
      chisq:      -999.0,$
      rms:        -999.0}
;     xaxis:    rr25axis,$
;     yfit:  rr25axis*0.0}

    for imonte = 0L, nmonte-1L do begin

       if (imonte eq 0) then xx = x else xx = x + randomn(seed,npts)*xerr
       if (imonte eq 0) then yy = y else yy = y + randomn(seed,npts)*yerr

       if (imonte eq 0) then begin
          if keyword_set(weightedfit) then begin
             thiscoeff = linfit(xx-xchar,yy,measure_errors=yerr,chisq=thischisq) 
          endif else begin
             thiscoeff = linfit(xx-xchar,yy,chisq=thischisq)
          endelse
          thiscoeff[0] = thiscoeff[0] - thiscoeff[1]*xchar
          thisyfit = poly(rr25axis,thiscoeff)
          yymodel = poly(xx,thiscoeff)
; store the best fit
;         res.yfit = thisyfit
          res.chisq = thischisq
          res.rms = djsig(yy-yymodel)
          res.int[0]   = thiscoeff[0]
          res.slope[0] = thiscoeff[1]
          res.ychar[0] = interpol(thisyfit,rr25axis,xchar)
          res.ycent[0] = interpol(thisyfit,rr25axis,xcent)
       endif 

       if (thiscoeff[1] gt 0.0) and (imonte eq 0L) then res.slope_flag = 1
       if (thiscoeff[1] gt 0.0) and (not keyword_set(noslopeconstrain)) then begin
          coeff = im_linefit(xx-xchar,yy,coeff_limits=[[0.0,0.0],[0.0,0.0]],$
            coeff_limited=[[0,0],[0,1]],coeff_guess=[8.5,-1.0]) ; slope must be negative
       endif else begin
          if keyword_set(weightedfit) then $
            coeff = linfit(xx-xchar,yy,measure_errors=yerr,chisq=chisq) else $
              coeff = linfit(xx-xchar,yy,chisq=chisq)
       endelse 

       coeff[0] = coeff[0] - coeff[1]*xchar
       yfit_monte = poly(rr25axis,coeff)

       if (imonte eq 0L) then begin
          int_monte = fltarr(nmonte)
          slope_monte = fltarr(nmonte)
          ychar_monte = fltarr(nmonte)
          ycent_monte = fltarr(nmonte)
       endif

       int_monte[imonte] = coeff[0]
       slope_monte[imonte] = coeff[1]
       ychar_monte[imonte] = interpol(yfit_monte,rr25axis,xchar)
       ycent_monte[imonte] = interpol(yfit_monte,rr25axis,xcent)

;      if (imonte eq 0L) then begin
;         plot, xx, yy, ps=4, xsty=3, ysty=3, sym=2, charsize=2
;         djs_oplot, rr25axis, thisyfit, color='red', thick=5
;      endif
;      oplot, rr25axis, yfit_monte
       
    endfor

    res.int[1]   = djsig(int_monte)
    res.slope[1] = djsig(slope_monte)
    res.ychar[1] = djsig(ychar_monte)
    res.ycent[1] = djsig(ycent_monte)

stop    
    
return, res
end

